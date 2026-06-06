#!/usr/bin/env Rscript
# h8_overcompression.R
# H8: The "paradox" — batches with jointly compressed biology AND noise get delta >> 1,
#     not delta < 1. ComBat then over-squashes them. Distinct from the flip (delta < 1) mode.
# Test:
#   1. Identify "compressed" batches: gap_raw_b < beta (weak biology) AND vl_b < vp (low noise)
#   2. Verify their delta.star >> 1 (the paradox)
#   3. Show their corrected class gap → beta (local biology erased, global β injected)
#   4. Contrast with "flip" batches (delta < 1): different profile

suppressMessages({ library(sva); library(matrixStats) })
setwd("/home/phr23/confounded_analysis/scripts/evaluations/book_chapter")
load("data/TB_real_data.RData")
OUT_DIR <- "outputs/diagnostics/hypothesis_tests"

ALL_STUDIES <- c("GSE37250_SA","USA","India","GSE37250_M","Africa","GSE39941_M")
bin <- function(l) as.integer(ifelse(l %in% c("1",1,"Active"), 1L, 0L))
train_studies <- function(n, test) ALL_STUDIES[ALL_STUDIES!=test][seq_len(n)]
log_safe <- function(m) { if(max(m,na.rm=TRUE)>100){mn<-min(m,na.rm=TRUE);if(mn<0)m<-m-mn;m<-log2(m+1)};m }

aprior <- function(gg) { m<-mean(gg); s2<-var(gg); (2*s2+m^2)/s2 }
bprior <- function(gg) { m<-mean(gg); s2<-var(gg); (m*s2+m^3)/s2 }
postvar_fn <- function(sum2, n, a, b) (0.5*sum2+b)/(n/2+a-1)

rows <- list()

for(nn in 2:5) for(tst in ALL_STUDIES) {
  refs <- train_studies(nn, tst)
  cg   <- Reduce(intersect, lapply(c(refs,tst), function(s) rownames(dat_lst[[s]])))
  dat  <- log_safe(do.call(cbind, lapply(refs, function(s) dat_lst[[s]][cg,])))
  bat  <- as.factor(do.call(c, lapply(refs, function(s) rep(s, ncol(dat_lst[[s]])))))
  lab  <- do.call(c, lapply(refs, function(s) bin(label_lst[[s]])))
  ubat <- levels(bat); N <- ncol(dat)

  # Global supervised LS fit
  des  <- cbind(model.matrix(~-1+bat), class=lab)
  B    <- solve(crossprod(des), crossprod(des, t(dat)))
  ce   <- B["class",]
  nbi  <- as.numeric(table(bat)); grand <- as.numeric(crossprod(nbi/N, B[seq_len(length(ubat)),]))
  vp   <- rowMeans((dat - t(des %*% B))^2)
  sm   <- outer(grand, rep(1,N)) + outer(ce, lab)
  Z    <- (dat - sm) / sqrt(pmax(vp, 1e-8))

  ref_sup <- suppressMessages(ComBat(dat, bat, mod=model.matrix(~lab)))
  w_sup   <- rowMeans(ref_sup[,lab==1,drop=FALSE]) - rowMeans(ref_sup[,lab==0,drop=FALSE])
  w_raw   <- rowMeans(dat[,lab==1,drop=FALSE])      - rowMeans(dat[,lab==0,drop=FALSE])

  for(k in seq_along(ubat)) {
    b  <- ubat[k]; i <- bat==b; n_b <- sum(i)
    n1 <- sum(i & lab==1); n0 <- sum(i & lab==0)
    if(n1<1 || n0<1) next

    gap_b  <- rowMeans(dat[,i & lab==1,drop=FALSE]) - rowMeans(dat[,i & lab==0,drop=FALSE])
    v1 <- if(n1>1) rowVars(dat[,i & lab==1,drop=FALSE]) else rep(0, nrow(dat))
    v0 <- if(n0>1) rowVars(dat[,i & lab==0,drop=FALSE]) else rep(0, nrow(dat))
    vl_b <- pmax((pmax(n1-1,1)*v1 + pmax(n0-1,1)*v0) / pmax(n1+n0-2,1), 1e-8)

    Zi   <- Z[,i,drop=FALSE]; ghat_b <- rowMeans(Zi)
    dr_b <- rowVars(Zi - ghat_b)
    if(n_b > 2 && var(dr_b) > 1e-10) {
      a_b <- aprior(dr_b); s_b <- bprior(dr_b)
      ds_b <- postvar_fn(dr_b*(n_b-1), n_b, a_b, s_b)
    } else { ds_b <- dr_b }

    # Corrected class gap for this batch (from supervised ComBat output)
    gap_sup_b <- rowMeans(ref_sup[,i & lab==1,drop=FALSE]) - rowMeans(ref_sup[,i & lab==0,drop=FALSE])

    # Classify each gene:
    # "compressed": gap_b < beta AND vl_b < vp   → expect delta >> 1
    # "flip":       delta_raw < 1                → expect delta < 1, sign inversion
    # "normal":     neither
    compressed <- gap_b < ce & vl_b < vp & abs(ce) > 0.01
    flip_zone  <- dr_b < 0.8  # conservative threshold for delta approaching flip zone

    # ── The paradox test: compressed batches should have large delta ──────────
    # Predicted delta² from the discussion formula: (mismatch + noise) / vp
    var_x <- (n1/n_b) * (n0/n_b)
    mismatch_b <- (gap_b - ce)^2 * var_x
    pred_delta2 <- (mismatch_b + vl_b) / pmax(vp, 1e-8)

    # For compressed genes: what is their actual delta.star?
    # The paradox predicts: although vl_b < vp, mismatch (biology−pooled)² blows up
    # so delta ends up >> 1.
    if(sum(compressed, na.rm=TRUE) > 20) {
      med_ds_compressed <- median(ds_b[compressed], na.rm=TRUE)
      med_dr_compressed <- median(dr_b[compressed], na.rm=TRUE)
      # What fraction of compressed genes have delta > 1? (paradox prediction: most)
      pct_ds_gt1_compressed <- 100*mean(ds_b[compressed] > 1, na.rm=TRUE)
      # What is the corrected class gap for these genes?
      med_gap_raw_c  <- median(gap_b[compressed],     na.rm=TRUE)
      med_gap_sup_c  <- median(gap_sup_b[compressed], na.rm=TRUE)
      med_beta_c     <- median(ce[compressed],         na.rm=TRUE)
      # Paradox metric: corrected gap for compressed genes → beta?
      gap_sup_vs_beta <- cor(gap_sup_b[compressed], ce[compressed], method="spearman", use="complete.obs")
      gap_raw_vs_beta <- cor(gap_b[compressed],     ce[compressed], method="spearman", use="complete.obs")
    } else {
      med_ds_compressed <- med_dr_compressed <- pct_ds_gt1_compressed <- NA
      med_gap_raw_c <- med_gap_sup_c <- med_beta_c <- NA
      gap_sup_vs_beta <- gap_raw_vs_beta <- NA
    }

    rows[[length(rows)+1]] <- data.frame(
      n=nn, test=tst, batch=b, n_b=n_b,
      # Batch-level profile
      pct_compressed = 100*mean(compressed, na.rm=TRUE),
      pct_flip_zone  = 100*mean(flip_zone,  na.rm=TRUE),
      med_gap_b      = median(gap_b,    na.rm=TRUE),
      med_ce         = median(ce,        na.rm=TRUE),
      med_vl_b       = median(vl_b,     na.rm=TRUE),
      med_vp         = median(vp,       na.rm=TRUE),
      med_dr_b       = median(dr_b,     na.rm=TRUE),
      med_ds_b       = median(ds_b,     na.rm=TRUE),
      # Paradox test (for compressed genes)
      pct_compressed_delta_gt1 = pct_ds_gt1_compressed,
      med_ds_compressed        = med_ds_compressed,
      med_gap_raw_compressed   = med_gap_raw_c,
      med_gap_sup_compressed   = med_gap_sup_c,
      med_beta_compressed      = med_beta_c,
      cor_gap_sup_vs_beta_compressed = gap_sup_vs_beta,
      cor_gap_raw_vs_beta_compressed = gap_raw_vs_beta,
      # Whole-batch gap alignment
      cos_sup_raw = {g1<-gap_sup_b;g2<-gap_b;sum(g1*g2)/(sqrt(sum(g1^2))*sqrt(sum(g2^2)))},
      cos_sup_beta= {g1<-gap_sup_b;b1<-ce;  sum(g1*b1)/(sqrt(sum(g1^2))*sqrt(sum(b1^2)))},
      stringsAsFactors=FALSE)
  }
  cat("done n=",nn," test=",tst,"\n",sep="")
}

res <- do.call(rbind, rows)
write.csv(res, file.path(OUT_DIR,"h8_overcompression.csv"), row.names=FALSE)

cat("\n=== H8: The over-compression paradox ===\n")
cat("Across all 24 scenarios x batches:\n")
cat("Median % genes 'compressed' (gap<beta AND vl<vp):",
    round(median(res$pct_compressed, na.rm=TRUE),1),"%\n")
cat("Median % genes in flip zone (delta<0.8):",
    round(median(res$pct_flip_zone, na.rm=TRUE),1),"%\n\n")
comp_batches <- res[!is.na(res$pct_compressed_delta_gt1),]
cat("Batches with ≥20 compressed genes:", nrow(comp_batches), "\n")
cat("Paradox: % of compressed genes with delta > 1:",
    round(median(comp_batches$pct_compressed_delta_gt1, na.rm=TRUE),1), "%\n")
cat("Compressed genes — corrected gap vs beta (Spearman r):",
    round(median(comp_batches$cor_gap_sup_vs_beta_compressed, na.rm=TRUE),3),
    " (should be high if biology erased -> global beta injected)\n")
cat("Compressed genes — raw gap vs beta (Spearman r):",
    round(median(comp_batches$cor_gap_raw_vs_beta_compressed, na.rm=TRUE),3), "\n")

cat("\n=== Breakdown: flip zone vs over-compression ===\n")
cat("Are flip-zone batches (low delta) different from compressed batches (expect high delta)?\n")
for(nn_v in 2:5) {
  sub <- res[res$n==nn_v,]
  cat(sprintf("n=%d: flip_zone median %.1f%%, compressed median %.1f%%\n",
      nn_v, median(sub$pct_flip_zone,na.rm=TRUE), median(sub$pct_compressed,na.rm=TRUE)))
}
cat("-> h8_overcompression.csv\n")
