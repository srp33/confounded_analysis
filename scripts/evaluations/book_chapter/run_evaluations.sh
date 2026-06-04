#!/bin/bash

#SBATCH --time 4:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=16G
#SBATCH -J "run_evaluations"
#SBATCH -o logs/evaluations_%A.log

BOOK_CHAPTER_DIR="$HOME/confounded_analysis/scripts/evaluations/book_chapter"
OUTPUT_DIR="$HOME/confounded_analysis/grp_batch_effects/outputs/book_chapter"

cd "$BOOK_CHAPTER_DIR"

mkdir -p "$OUTPUT_DIR/diagnostics/hypothesis_tests"

run_script() {
    local script="$1"
    echo ""
    echo "======================================================"
    echo "Running: $script"
    echo "Start: $(date)"
    echo "======================================================"
    pixi run Rscript "scripts/$script"
    local rc=$?
    echo "======================================================"
    echo "Finished: $script (exit code $rc) at $(date)"
    echo "======================================================"
    if [ $rc -ne 0 ]; then
        echo "ERROR: $script failed with exit code $rc. Aborting."
        exit $rc
    fi
}

# Step 1: master purity/alignment table (required by all downstream scripts)
run_script evaluate_mechanism.R

# Step 2: parallel-independent scripts (all read m1_pc1_alignment.csv)
run_script evaluate_vp_mechanism.R
run_script evaluate_mechanism_v2.R
run_script evaluate_targeted_compression.R
run_script evaluate_targeted_v2.R

# Step 3: hypothesis verdicts (reads all prior outputs)
run_script evaluate_hypotheses.R

echo ""
echo "All evaluation scripts completed successfully."
