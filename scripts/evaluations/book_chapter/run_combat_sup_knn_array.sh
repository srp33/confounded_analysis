#!/bin/bash
# run_combat_sup_knn_array.sh
# Re-run production combat_sup x knn across the 24 scenarios (n=2..5 x 6 test studies)
# as a SLURM array, after the step-2 leakage fix (mod=NULL).
#
# Usage:
#   sbatch run_combat_sup_knn_array.sh
# A dependent aggregation job (submitted separately with --dependency=afterok)
# merges the per-task CSVs and computes the mcc summary.

#SBATCH --array=0-23
#SBATCH --time 0:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=16G
#SBATCH -J "csup_knn"
#SBATCH -o logs/csup_knn_%A_%a.log

BOOK_CHAPTER_DIR="$HOME/confounded_analysis/scripts/evaluations/book_chapter"
cd "$BOOK_CHAPTER_DIR"
mkdir -p logs
OUTDIR="$BOOK_CHAPTER_DIR/outputs/diagnostics/combat_sup_knn_rerun"
mkdir -p "$OUTDIR"

STUDIES=(GSE37250_SA USA India GSE37250_M Africa GSE39941_M)
NVALS=(2 3 4 5)

i=${SLURM_ARRAY_TASK_ID}
n=${NVALS[$(( i / 6 ))]}
s=${STUDIES[$(( i % 6 ))]}
out="$OUTDIR/combat_sup_knn_n${n}_${s}.csv"

echo "Array task $i -> combat_sup x knn  n=$n test=$s"
echo "Start: $(date)"
pixi run Rscript scripts/classify_adjusters.R \
  --adjuster combat_sup --classifier knn \
  --num-datasets "$n" --test-study "$s" \
  -o "$out"
rc=$?
echo "Finished n=$n test=$s (exit $rc) at $(date)"
exit $rc
