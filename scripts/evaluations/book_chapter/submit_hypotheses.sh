#!/bin/bash
# submit_hypotheses.sh
# Submit the H1–H14 job array and a dependent aggregation job, then exit.
# Run this from the login node:  ./submit_hypotheses.sh

set -euo pipefail
cd "$(dirname "$0")"
mkdir -p logs

ARRAY_ID=$(sbatch --parsable run_hypotheses_array.sh)
echo "Submitted hypothesis array job: $ARRAY_ID (tasks 1-14)"

AGG_ID=$(sbatch --parsable \
  --dependency=afterok:"$ARRAY_ID" \
  --time=10:00 --ntasks=1 --cpus-per-task=1 --mem-per-cpu=8G \
  -J "hyp_aggregate" -o "logs/hyp_aggregate_%A.log" \
  --wrap "cd $PWD && pixi run Rscript scripts/aggregate_verdicts.R")
echo "Submitted aggregation job:    $AGG_ID (runs after array completes)"
echo ""
echo "Monitor with:  squeue -u $USER"
echo "Verdicts land in: outputs/diagnostics/hypothesis_tests/verdict_table.csv"
