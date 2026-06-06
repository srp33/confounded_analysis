#!/usr/bin/env Rscript
# evaluate_hypotheses.R
# Sequential driver for the H1–H14 evaluations. Sources the shared setup once,
# then runs each per-hypothesis script in scripts/hypotheses/, and aggregates
# the verdicts into verdict_table.csv.
#
# For a parallel run on the HPC use the SLURM job array instead:
#   ./submit_hypotheses.sh        (array of 14 + dependent aggregation)
# Both paths share scripts/hypotheses_common.R and scripts/hypotheses/H*.R,
# so the logic stays in one place.

source("scripts/hypotheses_common.R")

for (k in 1:14) {
  f <- sprintf("scripts/hypotheses/H%d.R", k)
  source(f)            # `sc` already exists, so the guard in H*.R is a no-op
}

source("scripts/aggregate_verdicts.R")
