#!/bin/bash
# aggregate_combat_sup_knn.sh
# Dependent job: merge the 24 combat_sup x knn per-scenario CSVs and summarize mcc.
#SBATCH --time 0:10:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
#SBATCH -J "csup_knn_agg"
#SBATCH -o logs/csup_knn_agg_%A.log

BOOK_CHAPTER_DIR="$HOME/confounded_analysis/scripts/evaluations/book_chapter"
cd "$BOOK_CHAPTER_DIR"
OUTDIR="$BOOK_CHAPTER_DIR/outputs/diagnostics/combat_sup_knn_rerun"

pixi run Rscript -e '
outdir <- file.path(Sys.getenv("BOOK_CHAPTER_DIR"), "outputs/diagnostics/combat_sup_knn_rerun")
files <- list.files(outdir, pattern="combat_sup_knn_.*\\.csv$", full.names=TRUE)
d <- do.call(rbind, lapply(files, function(f) tryCatch(read.csv(f), error=function(e) NULL)))
m <- d[tolower(d$metric)=="mcc",]
cat(sprintf("scenarios with mcc: %d / 24\n", nrow(m)))
cat("=== combat_sup x knn (post-fix, unsupervised step-2) ===\n")
cat(sprintf("  median MCC = %.3f\n", median(m$value, na.rm=TRUE)))
cat(sprintf("  mean   MCC = %.3f\n", mean(m$value, na.rm=TRUE)))
cat(sprintf("  below chance (<0) = %d / %d\n", sum(m$value<0, na.rm=TRUE), nrow(m)))
agg <- m[order(m$n_datasets, m$test_study), c("n_datasets","test_study","value")]
write.csv(agg, file.path(outdir, "AGG_mcc.csv"), row.names=FALSE)
print(agg, row.names=FALSE)
'
