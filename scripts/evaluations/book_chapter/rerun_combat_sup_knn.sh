#!/bin/bash
# Re-run production combat_sup x knn across the 24 scenarios (n=2..5 x 6 test studies)
# after the step-2 leakage fix (mod=NULL). Aggregates the mcc metric.
set -u
cd "$HOME/confounded_analysis/scripts/evaluations/book_chapter"
OUTDIR="/tmp/combat_sup_knn_rerun"
mkdir -p "$OUTDIR"
STUDIES=(GSE37250_SA USA India GSE37250_M Africa GSE39941_M)

for n in 2 3 4 5; do
  for s in "${STUDIES[@]}"; do
    out="$OUTDIR/combat_sup_knn_n${n}_${s}.csv"
    pixi run Rscript scripts/classify_adjusters.R \
      --adjuster combat_sup --classifier knn \
      --num-datasets "$n" --test-study "$s" \
      -o "$out" > "$OUTDIR/log_n${n}_${s}.txt" 2>&1
    echo "done n=$n test=$s -> exit $?"
  done
done

echo "=== Aggregating mcc across scenarios ==="
pixi run Rscript -e '
files <- list.files("/tmp/combat_sup_knn_rerun", pattern="combat_sup_knn_.*\\.csv$", full.names=TRUE)
d <- do.call(rbind, lapply(files, function(f) tryCatch(read.csv(f), error=function(e) NULL)))
m <- d[tolower(d$metric)=="mcc",]
cat(sprintf("scenarios with mcc: %d\n", nrow(m)))
cat(sprintf("combat_sup x knn (post-fix, unsupervised step-2):\n"))
cat(sprintf("  median MCC = %.3f\n", median(m$value, na.rm=TRUE)))
cat(sprintf("  mean   MCC = %.3f\n", mean(m$value, na.rm=TRUE)))
cat(sprintf("  n below chance (<0) = %d / %d\n", sum(m$value<0, na.rm=TRUE), nrow(m)))
write.csv(m[order(m$n_datasets, m$test_study), c("n_datasets","test_study","value")],
          "/tmp/combat_sup_knn_rerun/AGG_mcc.csv", row.names=FALSE)
'
echo "=== ALL DONE ==="
