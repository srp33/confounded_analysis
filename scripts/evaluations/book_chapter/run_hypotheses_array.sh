#!/bin/bash
# run_hypotheses_array.sh
# Evaluate H1–H14 in parallel as a SLURM job array, then aggregate the verdicts.
#
# Usage:
#   sbatch run_hypotheses_array.sh
#
# Each array task runs ONE hypothesis (scripts/hypotheses/H<id>.R), which sources
# the shared setup (scripts/hypotheses_common.R) and writes verdict_H<id>.csv.
# A dependent aggregation job merges them into verdict_table.csv once all finish.

#SBATCH --array=1-14
#SBATCH --time 1:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=16G
#SBATCH -J "hypotheses"
#SBATCH -o logs/hyp_%A_%a.log

BOOK_CHAPTER_DIR="$HOME/confounded_analysis/scripts/evaluations/book_chapter"
cd "$BOOK_CHAPTER_DIR"
mkdir -p logs

H="scripts/hypotheses/H${SLURM_ARRAY_TASK_ID}.R"
echo "======================================================"
echo "Array task ${SLURM_ARRAY_TASK_ID}: $H"
echo "Start: $(date)"
echo "======================================================"
pixi run Rscript "$H"
rc=$?
echo "======================================================"
echo "Finished $H (exit code $rc) at $(date)"
echo "======================================================"
exit $rc
