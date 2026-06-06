#!/usr/bin/env Rscript
# h6_h7_h9_delta_decomposition.R
#
# H6: Which pathway drives delta < 1?
#   "Mismatch Deficit" (small Dβ, outlier inflates denominator)
#   vs "Noise Deficit" (low native sigma²_true)
# H7: Option2 works because vl_b ≈ sigma²_true; ComBat's delta_raw mixes in Dβ²·Var(X).
#   Validation: delta_raw · vp ≈ vl + (gap_raw - beta)² · Var(X)  [Var(X)=p·(1-p)]
# H9: The decomposition formula predicts which genes flip,
#   better than just delta_raw < threshold.

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

rows_batch <- list()   # per-batch summaries
rows_h9    <- list()   # per-scenario flip prediction quality

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

  # Per-batch quantities
  for(k in seq_along(ubat)) {
    b  <- ubat[k]; i  <- bat == b; n_b <- sum(i)
    n1 <- sum(i & lab==1); n0 <- sum(i & lab==0)
    if(n1 < 1 || n0 < 1) next
    p_b <- n1 / n_b  # class proportion (Var(X) = p*(1-p))

    # Raw class gap in this batch
    gap_b <- rowMeans(dat[,i & lab==1,drop=FALSE]) - rowMeans(dat[,i & lab==0,drop=FALSE])

    # Noise component: within-class pooled variance (sigma²_true proxy)
    v1  <- if(n1>1) rowVars(dat[,i & lab==1,drop=FALSE]) else rep(0, nrow(dat))
    v0  <- if(n0>1) rowVars(dat[,i & lab==0,drop=FALSE]) else rep(0, nrow(dat))
    vl_b <- pmax((pmax(n1-1,1)*v1 + pmax(n0-1,1)*v0) / pmax(n1+n0-2,1), 1e-8)

    # Mismatch component: (gap_b - beta)² · Var(X)
    var_x <- p_b * (1 - p_b)  # scalar per batch
    mismatch_b <- (gap_b - ce)^2 * var_x  # per gene

    # ComBat delta_raw = var(Z - gamma_hat | batch)
    Zi     <- Z[,i,drop=FALSE]; ghat_b <- rowMeans(Zi)
    dr_b   <- rowVars(Zi - ghat_b)

    # ── H7 validation: delta_raw · vp ≈ vl_b + mismatch_b ──────────────────
    # In the raw expression scale: ComBat's local variance ≈ noise + mismatch
    # delta_raw * vp is ComBat's estimate of local variance in raw space
    combat_local_var <- dr_b * vp
    predicted_var    <- vl_b + mismatch_b

    cor_h7  <- cor(combat_local_var, predicted_var, method="spearman", use="complete.obs")
    # R² of linear fit through origin: predicted_var ~ combat_local_var
    lm_coef <- sum(predicted_var * combat_local_var, na.rm=TRUE) / sum(combat_local_var^2, na.rm=TRUE)
    resid_ss <- sum((predicted_var - lm_coef * combat_local_var)^2, na.rm=TRUE)
    tot_ss   <- sum((predicted_var - mean(predicted_var, na.rm=TRUE))^2, na.rm=TRUE)
    r2_h7    <- 1 - resid_ss/tot_ss

    # ── H6: pathway classification for genes where delta_raw < 1 ────────────
    below1 <- dr_b < 1
    if(sum(below1, na.rm=TRUE) > 0) {
      noise_frac     <- median(vl_b[below1]     / pmax(vp[below1], 1e-8), na.rm=TRUE)
      mismatch_frac  <- median(mismatch_b[below1] / pmax(vp[below1], 1e-8), na.rm=TRUE)
    } else {
      noise_frac <- mismatch_frac <- NA
    }
    # Also: among delta_raw < 1 genes, is vl_b < vp (Noise Deficit) or mismatch near 0?
    noise_deficit_pct    <- if(sum(below1,na.rm=TRUE)>0)
      100*mean(vl_b[below1] < vp[below1], na.rm=TRUE) else NA
    mismatch_deficit_pct <- if(sum(below1,na.rm=TRUE)>0)
      100*mean(mismatch_b[below1] < 0.1*vp[below1], na.rm=TRUE) else NA

    # ── EB-shrunk delta.star for this batch ──────────────────────────────────
    if(n_b > 2 && var(dr_b) > 1e-10) {
      a_b <- aprior(dr_b); s_b <- bprior(dr_b)
      ds_b <- postvar_fn(dr_b*(n_b-1), n_b, a_b, s_b)
    } else { ds_b <- dr_b }

    rows_batch[[length(rows_batch)+1]] <- data.frame(
      n=nn, test=tst, batch=b, n_b=n_b, n1=n1, n0=n0,
      # H6: pathway breakdown
      pct_below1_draw     = 100*mean(dr_b < 1, na.rm=TRUE),
      pct_below1_dstar    = 100*mean(ds_b < 1, na.rm=TRUE),
      noise_deficit_pct   = noise_deficit_pct,    # among below1: % with vl < vp
      mismatch_deficit_pct= mismatch_deficit_pct, # among below1: % with mismatch ≈ 0
      noise_frac_median   = noise_frac,
      mismatch_frac_median= mismatch_frac,
      # H7: formula validation
      cor_h7  = cor_h7,
      r2_h7   = r2_h7,
      lm_coef = lm_coef,
      # Medians for context
      med_dr   = median(dr_b), med_ds = median(ds_b),
      med_vl   = median(vl_b), med_vp = median(vp),
      med_mismatch = median(mismatch_b),
      stringsAsFactors=FALSE)
  }

  # ── H9: per-scenario flip prediction quality ─────────────────────────────
  # For each gene+batch, compare:
  #   "analytic prediction": delta_raw < flip_thr (uses observed delta_raw)
  #   "formula prediction":  predicted_var/vp < flip_thr  (uses vl + mismatch)
  #   "empirical flip":      sign(corrected_gap) ≠ sign(raw_gap)
  ref_sup <- suppressMessages(ComBat(dat, bat, mod=model.matrix(~lab)))
  w_raw <- rowMeans(dat[,lab==1,drop=FALSE])    - rowMeans(dat[,lab==0,drop=FALSE])
  w_sup <- rowMeans(ref_sup[,lab==1,drop=FALSE]) - rowMeans(ref_sup[,lab==0,drop=FALSE])
  empirical_flip <- sign(w_sup) != sign(w_raw) & abs(w_raw) > 0.01

  # Per-batch: compute flip_thr and predicted_delta
  for(k in seq_along(ubat)) {
    b <- ubat[k]; i <- bat==b; n_b <- sum(i)
    n1 <- sum(i & lab==1); n0 <- sum(i & lab==0)
    if(n1<1 || n0<1) next
    p_b <- n1/n_b
    gap_b <- rowMeans(dat[,i & lab==1,drop=FALSE]) - rowMeans(dat[,i & lab==0,drop=FALSE])
    v1  <- if(n1>1) rowVars(dat[,i & lab==1,drop=FALSE]) else rep(0, nrow(dat))
    v0  <- if(n0>1) rowVars(dat[,i & lab==0,drop=FALSE]) else rep(0, nrow(dat))
    vl_b <- pmax((pmax(n1-1,1)*v1 + pmax(n0-1,1)*v0) / pmax(n1+n0-2,1), 1e-8)
    mismatch_b <- (gap_b - ce)^2 * p_b * (1-p_b)
    pred_var_b <- (vl_b + mismatch_b) / pmax(vp, 1e-8)  # predicted delta²

    Zi <- Z[,i,drop=FALSE]; ghat_b <- rowMeans(Zi)
    dr_b <- rowVars(Zi - ghat_b)

    # Flip threshold per gene: only defined where gap_raw < beta and beta > 0
    valid <- ce > 0 & gap_b < ce & abs(ce) > 0.01
    if(sum(valid) < 50) next
    ft <- 1 - gap_b[valid]/ce[valid]  # ComBat flip threshold (delta < ft → flip)

    # Predicted flip (using formula): pred_var < ft (= flip_thr for ComBat)
    flip_by_draw   <- dr_b[valid] < ft
    flip_by_formula<- pred_var_b[valid] < ft

    rows_h9[[length(rows_h9)+1]] <- data.frame(
      n=nn, test=tst, batch=b,
      # How well does delta_raw predict empirical flips?
      cor_draw_flip    = cor(as.integer(dr_b[valid] < ft), empirical_flip[valid], method="pearson"),
      # How well does the formula predict delta_raw < ft?
      agree_formula_draw= mean(flip_by_formula == flip_by_draw, na.rm=TRUE),
      # Does the formula's predicted delta correlate with observed delta?
      cor_formula_draw  = cor(pred_var_b[valid], dr_b[valid], method="spearman"),
      n_valid=sum(valid), n_flips_empirical=sum(empirical_flip[valid],na.rm=TRUE),
      n_pred_draw=sum(flip_by_draw,na.rm=TRUE),
      n_pred_formula=sum(flip_by_formula,na.rm=TRUE),
      stringsAsFactors=FALSE)
  }

  cat("done n=",nn," test=",tst,"\n",sep="")
}

rb <- do.call(rbind, rows_batch)
r9 <- do.call(rbind, rows_h9)
write.csv(rb, file.path(OUT_DIR,"h6_h7_delta_decomposition.csv"), row.names=FALSE)
write.csv(r9, file.path(OUT_DIR,"h9_flip_prediction.csv"), row.names=FALSE)

cat("\n=== H7: formula validation delta_raw·vp ≈ vl_b + mismatch_b ===\n")
cat("Median Spearman r:", round(median(rb$cor_h7, na.rm=TRUE),3), "\n")
cat("Median R²:        ", round(median(rb$r2_h7,  na.rm=TRUE),3), "\n")
cat("Median lm coef (should be ~1):", round(median(rb$lm_coef, na.rm=TRUE),3), "\n\n")

cat("=== H6: Pathway to delta < 1 ===\n")
below1_batches <- rb[rb$pct_below1_draw > 0 & !is.na(rb$pct_below1_draw),]
cat("Batches with any delta_raw < 1:", nrow(below1_batches), "of", nrow(rb), "\n")
cat("Among those: median noise_deficit_pct (vl<vp among below1 genes):",
    round(median(below1_batches$noise_deficit_pct, na.rm=TRUE),1), "%\n")
cat("Among those: median mismatch_deficit_pct (mismatch≈0 among below1 genes):",
    round(median(below1_batches$mismatch_deficit_pct, na.rm=TRUE),1), "%\n")
cat("-> dominant pathway is whichever % is higher\n\n")

cat("=== H9: Formula predicts flips ===\n")
cat("Median agreement (formula vs delta_raw prediction):",
    round(median(r9$agree_formula_draw, na.rm=TRUE),3), "\n")
cat("Median Spearman cor (formula delta² vs actual delta_raw):",
    round(median(r9$cor_formula_draw, na.rm=TRUE),3), "\n")
cat("-> h6_h7_delta_decomposition.csv, h9_flip_prediction.csv\n")
