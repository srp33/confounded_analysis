#!/usr/bin/env python3
"""
Check for failed Snakemake jobs by comparing log files to expected outputs.
Also extracts error messages from failed logs.
"""

import sys
from pathlib import Path
import re

def load_config(config_file):
    """Load the Snakemake config - simple YAML parser"""
    config = {}
    try:
        with open(config_file) as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith('#') and ':' in line:
                    key, value = line.split(':', 1)
                    key = key.strip()
                    value = value.strip().strip('"\'')
                    config[key] = value
    except Exception as e:
        print(f"Warning: Could not parse config: {e}")
    return config

def extract_error_from_log(log_file, num_lines=50):
    """Extract error message from a log file, looking for actual errors"""
    try:
        with open(log_file) as f:
            lines = f.readlines()
        
        # Look for error patterns in the last N lines
        error_patterns = [
            r'Error',
            r'ERROR',
            r'error',
            r'Execution halted',
            r'Traceback',
            r'Exception',
            r'FAILED',
            r'failed'
        ]
        
        # Get last num_lines
        tail_lines = lines[-num_lines:] if len(lines) > num_lines else lines
        
        # Find the first line with an error pattern
        error_start_idx = None
        for i, line in enumerate(tail_lines):
            if any(re.search(pattern, line) for pattern in error_patterns):
                error_start_idx = i
                break
        
        # If we found an error, show from that point to the end
        if error_start_idx is not None:
            error_lines = tail_lines[error_start_idx:]
            # Remove empty lines at the end
            while error_lines and not error_lines[-1].strip():
                error_lines.pop()
            # Take up to 15 lines from the error point
            error_lines = error_lines[:15]
            return "".join(error_lines).rstrip()
        
        # If no error pattern found, just show the last 10 non-empty lines
        non_empty = [l for l in tail_lines if l.strip()]
        return "".join(non_empty[-10:]).rstrip()
        
    except Exception as e:
        return f"Could not read log: {e}"

def check_failures(output_folder):
    """Check for jobs that have logs but no outputs"""
    failures = []
    
    checks = [
        {
            "name": "classify_adjusters",
            "log_dir": Path(output_folder) / "logs" / "classify_adjusters",
            "result_dir": Path(output_folder) / "results" / "adjusters" / "individual"
        },
        {
            "name": "classify_batch_effects",
            "log_dir": Path(output_folder) / "logs" / "classify_batch_effects",
            "result_dir": Path(output_folder) / "results" / "batch_effects" / "individual"
        }
    ]
    
    for check in checks:
        if check["log_dir"].exists():
            for log_file in check["log_dir"].glob("*.log"):
                expected_output = check["result_dir"] / f"{log_file.stem}.csv"
                if not expected_output.exists():
                    failures.append({
                        "rule": check["name"],
                        "log": str(log_file),
                        "expected_output": str(expected_output),
                        "error": extract_error_from_log(log_file)
                    })
    
    return failures

def main():
    # Allow override via command line
    if len(sys.argv) > 1:
        output_folder = sys.argv[1]
    else:
        # Load config
        script_dir = Path(__file__).parent
        config_file = script_dir / "config.yaml"
        
        if not config_file.exists():
            print(f"Error: Config file not found at {config_file}")
            sys.exit(1)
        
        config = load_config(config_file)
        output_folder = config.get("output_folder", "outputs/book_chapter")
    
    output_folder = Path(output_folder)
    
    if not output_folder.exists():
        print(f"Error: Output folder does not exist: {output_folder}")
        print(f"\nTip: Run this script from the repo root, or provide the output folder path:")
        print(f"  python3 scripts/evaluations/book_chapter/check_failures.py <output_folder>")
        sys.exit(1)
    
    print(f"Checking for failures in: {output_folder}")
    print("=" * 80)
    
    failures = check_failures(output_folder)
    
    if not failures:
        print("\n✓ No failures found! All jobs with logs have corresponding outputs.\n")
        return 0
    
    print(f"\n⚠️  Found {len(failures)} failed jobs:\n")
    
    # Group by error type - extract key error message
    def get_error_signature(error_text):
        """Extract a signature from error text for grouping"""
        # Look for common error patterns
        patterns = [
            (r'Error: (.+?)(?:\n|$)', 'Error'),
            (r'ERROR\] Error message: (.+?)(?:\n|$)', 'Error'),
            (r'Error in (.+?) :', 'Error in'),
            (r'Execution halted', 'Execution halted'),
            (r"Couldn't find a sufficient Python binary", 'Python not found'),
            (r'Mean parameter must be a positive number', 'Invalid mean parameter'),
            (r"invalid 'row.names' length", 'Invalid row.names'),
            (r'Memory usage:', 'Memory/Unknown error'),
        ]
        
        for pattern, label in patterns:
            match = re.search(pattern, error_text)
            if match:
                if match.groups():
                    return f"{label}: {match.group(1)[:60]}"
                return label
        
        # If no pattern matches, use first line with "error" (case insensitive)
        for line in error_text.split('\n'):
            if re.search(r'error|Error|ERROR', line):
                return line[:80]
        
        # Last resort: first 80 chars
        return error_text[:80]
    
    error_groups = {}
    for fail in failures:
        error_key = get_error_signature(fail["error"])
        if error_key not in error_groups:
            error_groups[error_key] = []
        error_groups[error_key].append(fail)
    
    # Sort by number of failures (most common first)
    error_groups = dict(sorted(error_groups.items(), key=lambda x: len(x[1]), reverse=True))
    
    print(f"Grouped into {len(error_groups)} unique error type(s):\n")
    
    for i, (error_key, group) in enumerate(error_groups.items(), 1):
        print(f"\n{'='*80}")
        print(f"Error Type {i}: {len(group)} jobs failed with similar error")
        print(f"{'='*80}")
        print(f"\nSample error from: {group[0]['log']}")
        print(f"\n{group[0]['error']}")
        print(f"\nAffected jobs ({len(group)} total):")
        for fail in group[:10]:  # Show first 10
            print(f"  - {Path(fail['log']).name}")
        if len(group) > 10:
            print(f"  ... and {len(group) - 10} more")
    
    print(f"\n{'='*80}")
    print(f"\nSummary: {len(failures)} total failures")
    print(f"See individual log files for full details.")
    print(f"{'='*80}\n")
    
    return 1

if __name__ == "__main__":
    sys.exit(main())
