#!/usr/bin/env Rscript
# aggregate_verdicts.R
# Merge the per-hypothesis verdict CSVs (verdict_H*.csv) written by the parallel
# hypothesis scripts into a single ordered verdict_table.csv.

OUT_DIR <- "outputs/diagnostics/hypothesis_tests"

files <- list.files(OUT_DIR, pattern = "^verdict_H[0-9]+\\.csv$", full.names = TRUE)
if (length(files) == 0) stop("No verdict_H*.csv files found in ", OUT_DIR)

VT <- do.call(rbind, lapply(files, function(f)
  read.csv(f, stringsAsFactors = FALSE)))

# Order by the numeric part of the H tag (H1, H2, ... H14)
hyp_num <- as.integer(sub(".*?([0-9]+).*", "\\1", VT$H))
VT <- VT[order(hyp_num), ]

cat("\n\n═══════════════════════════════════════════════════════════════════\n")
cat("  SUMMARY\n")
cat("═══════════════════════════════════════════════════════════════════\n")
print(VT, row.names = FALSE)

expected <- paste0("H", 1:14)
got      <- sub(" .*", "", VT$H)
missing  <- setdiff(expected, got)
if (length(missing) > 0)
  cat(sprintf("\n  WARNING: missing verdicts for %s\n", paste(missing, collapse = ", ")))

write.csv(VT, file.path(OUT_DIR, "verdict_table.csv"), row.names = FALSE)
cat(sprintf("\nVerdict table written to %s/verdict_table.csv (%d/14 hypotheses)\n",
            OUT_DIR, nrow(VT)))
