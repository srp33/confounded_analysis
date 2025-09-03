# Cache utilities - HashCache and DataFrameCache moved from scripts/metrics/util.py

import os
import hashlib
import json
import contextlib
import time
import pandas as pd
from pathlib import Path

class HashCache:
    """
    Manages file hashes to avoid re-computing results for unchanged data.
    Includes generalized logic to manage and "top up" repeated runs.
    """
    def __init__(self, hash_dir, cache_filename, write_over=False, debug=False):
        self.write_over = write_over
        self.debug = debug
        self.cache_path = Path(hash_dir) / cache_filename
        self.lock_path = Path(hash_dir) / f"{cache_filename}.lock"
        self.hashes = self._load_hashes()

    def _acquire_lock(self):
        """Acquires a file-based lock, waiting if necessary."""
        while True:
            try:
                os.mkdir(self.lock_path)
                return
            except FileExistsError:
                time.sleep(0.1)

    def _release_lock(self):
        """Releases the file-based lock."""
        try:
            os.rmdir(self.lock_path)
        except (FileNotFoundError, OSError) as e:
            if self.debug: print(f"DEBUG (PID {os.getpid()}): Error releasing lock: {e}")

    def _load_hashes(self):
        """Loads the hash dictionary from a JSON file, using a lock."""
        if self.write_over:
            if self.debug: print("DEBUG: Write-over enabled. Starting with an empty cache.")
            return {}
        
        self._acquire_lock()
        try:
            if not self.cache_path.exists():
                return {}
            with open(self.cache_path, 'r') as f:
                return json.load(f)
        except (FileNotFoundError, json.JSONDecodeError):
            if self.debug: print(f"DEBUG (PID {os.getpid()}): Hash file not found or invalid. Starting with empty cache.")
            return {}
        finally:
            self._release_lock()

    def _save_hashes(self):
        """Saves the current hash dictionary to a JSON file, using a lock."""
        self._acquire_lock()
        try:
            current_hashes_on_disk = {}
            if self.cache_path.exists():
                with open(self.cache_path, 'r') as f:
                    try:
                        current_hashes_on_disk = json.load(f)
                    except json.JSONDecodeError:
                        pass
            
            current_hashes_on_disk.update(self.hashes)
            self.hashes = current_hashes_on_disk

            self.cache_path.parent.mkdir(parents=True, exist_ok=True)
            with open(self.cache_path, 'w') as f:
                json.dump(self.hashes, f, indent=4)
            if self.debug:
                print(f"DEBUG (PID {os.getpid()}): Saved {len(self.hashes)} hashes to {self.cache_path}")
        finally:
            self._release_lock()

    @staticmethod
    def _calculate_hash(file_path):
        """Calculates the MD5 hash of a single file."""
        try:
            with open(file_path, 'rb') as f:
                file_hash = hashlib.md5()
                while chunk := f.read(8192):
                    file_hash.update(chunk)
            return file_hash.hexdigest()
        except FileNotFoundError:
            return None

    @staticmethod
    def _calculate_combined_hash(file_paths):
        """Calculates a single hash for a list of files."""
        combined_hasher = hashlib.md5()
        for file_path in sorted(file_paths):
            file_hash = HashCache._calculate_hash(file_path)
            if file_hash:
                combined_hasher.update(file_hash.encode())
        return combined_hasher.hexdigest()

    def _build_mask(self, df, identifier_dict):
        """Helper to build a pandas boolean mask from a run identifier dictionary."""
        # Start with a mask that is all True
        mask = pd.Series(True, index=df.index)
        if not identifier_dict:
            return ~mask # Return all False if no identifier
        # Iteratively apply filters for each key-value pair in the identifier
        for col, value in identifier_dict.items():
            if col not in df.columns:
                # If a key column is missing, this run can't exist in the file.
                raise KeyError(f"Identifier column '{col}' not found in DataFrame.")
            if isinstance(value, list):
                mask &= df[col].isin(value)
            else:
                mask &= df[col] == value
        return mask

    def _count_existing_repeats(self, output_file, run_identifier, count_config):
        """
        Counts completed repeats by checking results in the output file. A repeat
        is considered complete only if all its constituent metrics are present.
        """
        # Return early if required configuration is missing.
        if not all([output_file, run_identifier, count_config]):
            return 0
        try:
            df = pd.read_csv(output_file)
            
            # Extract configuration values to fail fast on bad config and for clarity.
            pg_col = count_config['primary_grouping_col']
            pg_val = count_config['primary_grouping_val']
            count_col = count_config['count_col']
            
            # Build the base filter mask for the specific experimental run.
            run_mask = self._build_mask(df, run_identifier)
            
            # Filter and count 
            value_counts = df[run_mask & (df[pg_col] == pg_val)][count_col].value_counts()
            
            # If value_counts is empty, no repeats exist. Otherwise, the number of
            # completed repeats is the minimum count across all tracked items.
            # Using .min() correctly identifies the number of *fully* completed repeats.
            return value_counts.min() if not value_counts.empty else 0
            
        except (FileNotFoundError, pd.errors.EmptyDataError, KeyError) as e:
            # Catch file I/O errors, empty files, or missing columns/config keys.
            if self.debug:
                print(f"DEBUG: Could not count repeats for {run_identifier}. Reason: {e}")
            return 0

    def _clear_existing_repeats(self, output_file, run_identifier):
        """Removes results for a specific run from the output file."""
        if not output_file or not run_identifier:
            return
        try:
            df = pd.read_csv(output_file)
            
            # Identify rows to remove based on the run identifier
            mask_to_remove = self._build_mask(df, run_identifier)
            
            # If no rows match the removal criteria, do nothing
            if not mask_to_remove.any():
                return

            print(f"Clearing stale results for run matching {run_identifier} from {output_file}")
            # Keep all rows *not* in the removal mask
            df[~mask_to_remove].to_csv(output_file, index=False)
        except (FileNotFoundError, pd.errors.EmptyDataError, KeyError) as e:
            if self.debug:
                print(f"DEBUG: Could not clear repeats for {run_identifier}. Reason: {e}")
            pass

    @contextlib.contextmanager
    def check(self, key, input_path):
        """
        Backwards-compatible check for simple skip/run logic.
        Yields a boolean: True to skip, False to run.
        """
        input_paths = [input_path] if isinstance(input_path, (str, Path)) else input_path

        # Pass None for repeat management parameters as they are not used here.
        with self.check_and_manage_repeats(
            key=key, 
            input_paths=input_paths, 
            n_repeats_requested=1, 
            output_file=None,
            run_identifier=None,
            count_config=None
        ) as (action, n_to_run, n_existing):
            should_skip = (action == "SKIP")
            yield should_skip

    @contextlib.contextmanager
    def check_and_manage_repeats(self, key, input_paths, n_repeats_requested, output_file, run_identifier, count_config):
        """
        Context manager to check hashes, manage partial runs, and clear stale results.
        Yields a status, number of repeats to run, and number of existing repeats.
        """
        current_hash = self._calculate_combined_hash(input_paths)
        previous_hash = self.hashes.get(key)
        hashes_match = (current_hash == previous_hash)
        
        action = "SKIP"
        n_to_run = 0
        n_existing = 0

        if self.write_over:
            print(f"Forcing re-run for '{key}'. Clearing old results.")
            self._clear_existing_repeats(output_file, run_identifier)
            action = "RUN_FULL"
            n_to_run = n_repeats_requested
        elif not hashes_match:
            print(f"Input data changed for '{key}'. Clearing old results.")
            self._clear_existing_repeats(output_file, run_identifier)
            action = "RUN_FULL"
            n_to_run = n_repeats_requested
        elif not output_file:
            # If no output file is specified, we can't manage repeats.
            # Since the hash check passed, we can assume a run was performed, and we can skip.
            action = "SKIP"
        else:
            n_existing = self._count_existing_repeats(output_file, run_identifier, count_config)
            if n_repeats_requested > n_existing:
                n_to_run = n_repeats_requested - n_existing
                print(f"Cache hit for '{key}'. Existing: {n_existing}. Topping up with {n_to_run} more repeats.")
                action = "RUN_PARTIAL"
            else:
                print(f"Skipping '{key}'. Found {n_existing} repeats, meeting requirement of {n_repeats_requested}.")
                action = "SKIP"
        
        try:
            yield action, n_to_run, n_existing
            if action == "RUN_FULL":
                self.hashes[key] = current_hash
        except Exception as e:
            print(f"ERROR during processing for key '{key}'. Hash will NOT be updated. Error: {e}")
            raise e
        finally:
            self._save_hashes()


class DataFrameCache(object):
    def __init__(self, folder:Path=None):
        self.dataframes = {} # {path: dataframe}
        self.folder = folder

    def __get_filepath__(self, file_name):
        if self.folder is None:
            raise ValueError("folder must be set to use file_name.")
        return os.path.join(self.folder, file_name)

    def __validate_file_path__(self, file_path, file_name):
        if file_path is None and file_name is None:
            raise ValueError("Either file_path or file_name must be provided.")
        if file_name is not None:
            file_path = self.__get_filepath__(file_name)
        return file_path

    def set_dataframe(self, df, file_path=None, file_name=None):
        file_path = self.__validate_file_path__(file_path, file_name)
        
        if not isinstance(df, pd.DataFrame):
            raise ValueError("df must be a pandas DataFrame.")
        
        self.dataframes[file_path] = df

    def get_dataframe(self, file_path=None, file_name=None):
        if file_name in self.dataframes:
            return self.dataframes[file_name]
            
        file_path = self.__validate_file_path__(file_path, file_name)
            
        if file_path in self.dataframes:
            return self.dataframes[file_path]
        else:
            try:
                df = pd.read_csv(file_path)
                self.dataframes[file_path] = df
                return df
            except FileNotFoundError:
                print(f"Process {os.getpid()}: ERROR - File not found: {file_path}", flush=True)
                return pd.DataFrame()
            except Exception as e:
                print(f"Process {os.getpid()}: ERROR - Could not read file {file_path} due to: {e}", flush=True)
                return pd.DataFrame()