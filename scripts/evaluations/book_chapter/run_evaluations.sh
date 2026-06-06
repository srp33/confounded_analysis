#!/bin/bash

#SBATCH --time 4:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=8G
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
        echo "ERROR: $script failed with exit code $rc."
        return $rc
    fi
}

run_script_bg() {
    local script="$1"
    local logfile="logs/eval_${script%.R}_$$.log"
    echo "Launching background: $script -> $logfile"
    (
        echo "======================================================"
        echo "Running: $script"
        echo "Start: $(date)"
        echo "======================================================"
        pixi run Rscript "scripts/$script"
        rc=$?
        echo "======================================================"
        echo "Finished: $script (exit code $rc) at $(date)"
        echo "======================================================"
        exit $rc
    ) > "$logfile" 2>&1 &
    echo $!
}

# Step 1: master purity/alignment table (required by all downstream scripts)
run_script evaluate_mechanism.R || exit $?

# Step 2: parallel-independent scripts (all read m1_pc1_alignment.csv)
pids=()
scripts_bg=(evaluate_vp_mechanism.R evaluate_mechanism_v2.R evaluate_targeted_compression.R evaluate_targeted_v2.R)
for s in "${scripts_bg[@]}"; do
    pid=$(run_script_bg "$s")
    pids+=("$pid:$s")
done

# Wait for all background scripts and check exit codes
failed=0
for entry in "${pids[@]}"; do
    pid="${entry%%:*}"
    script="${entry##*:}"
    wait "$pid"
    rc=$?
    if [ $rc -ne 0 ]; then
        echo "ERROR: $script failed with exit code $rc."
        failed=1
    else
        echo "OK: $script completed successfully."
    fi
done

if [ $failed -ne 0 ]; then
    echo "One or more parallel scripts failed. Aborting."
    exit 1
fi

# Merge background logs into main log
for s in "${scripts_bg[@]}"; do
    logfile="logs/eval_${s%.R}_$$.log"
    echo ""
    echo "=== Output from $s ==="
    cat "$logfile"
    rm -f "$logfile"
done

# Step 3: hypothesis verdicts (reads all prior outputs)
run_script evaluate_hypotheses.R || exit $?

echo ""
echo "All evaluation scripts completed successfully."
