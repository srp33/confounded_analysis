#!/usr/bin/env python3
"""
Run complete POSSE diagnostic analysis including timidity plots.
This script orchestrates the full diagnostic workflow.
"""

import subprocess
import sys
import os
from pathlib import Path

def run_command(cmd, description):
    """Run a command and handle errors."""
    print(f"\n=== {description} ===")
    print(f"Running: {' '.join(cmd)}")
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode != 0:
        print(f"ERROR: {description} failed")
        print(f"STDERR: {result.stderr}")
        return False
    else:
        print(f"SUCCESS: {description} completed")
        if result.stdout.strip():
            print(f"OUTPUT: {result.stdout}")
        return True

def main():
    """Run the complete diagnostic workflow."""
    
    # Define paths
    base_dir = Path("outputs")
    diagnostics_dir = base_dir / "diagnostics" / "posse"
    output_dir = base_dir
    
    print("=== POSSE v4.1 Complete Diagnostic Analysis ===")
    print(f"Base directory: {base_dir}")
    print(f"Diagnostics directory: {diagnostics_dir}")
    
    # Step 1: Aggregate POSSE diagnostics
    cmd1 = [
        "python", "scripts/aggregate_posse_diagnostics.py",
        "--input-dir", str(diagnostics_dir),
        "--output-summary", str(output_dir / "posse_diagnostic_summary.csv"),
        "--output-report", str(output_dir / "posse_diagnostic_report.md")
    ]
    
    if not run_command(cmd1, "POSSE Diagnostic Aggregation"):
        return 1
    
    # Step 2: Generate timidity analysis plots
    cmd2 = [
        "python", "scripts/generate_timidity_plots.py",
        "--output-dir", str(output_dir / "timidity_analysis")
    ]
    
    if not run_command(cmd2, "Timidity Analysis Generation"):
        return 1
    
    print("\n=== DIAGNOSTIC ANALYSIS COMPLETE ===")
    print(f"Results available in: {output_dir}")
    print("\nGenerated files:")
    print(f"  - {output_dir}/posse_diagnostic_summary.csv")
    print(f"  - {output_dir}/posse_diagnostic_report.md")
    print(f"  - {output_dir}/timidity_analysis/posse_vs_combat_timidity_analysis.png")
    print(f"  - {output_dir}/timidity_analysis/timidity_analysis_report.md")
    print(f"  - {output_dir}/timidity_analysis/correction_comparison_data.npz")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())