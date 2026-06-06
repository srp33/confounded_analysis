#!/usr/bin/env Rscript
# h3_h4_cohens_d_gap.R
# H3: Option2 preserves Cohen's d within each batch; mean.only does not.
# H4: Option2 and mean.only produce near-identical pooled class gap vectors (cosine sim).

suppressMessages({ library(sva); library(matrixStats) })
setwd("/home/phr23/confounded_analysis/scripts/evaluations/book_chapter")
load("data/TB_real_data.RData")
OUT_DIR <- "outputs/diagnostics/hypothesis_tests"

ALL_STUDIES <- c("GSE37250_SA","USA","India","GSE37250_M","Africa","GSE39941_M")
bin <- function(l) as.integer(ifelse(l %in% c("1",1,"Active"), 1L, 0L))
train_studies <- function(n, test) ALL_STUDIES[ALL_STUDIES!=test][seq_len(n)]
log_safe <- function(m) { if(max(m,na.rm=TRUE)>100){ mn<-min(m,na.rm=TRUE); if(mn<0) m<-m-mn; m<-log2(m+1) }; m }
cosim <- function(a,b) sum(a*b)/(sqrt(sum(a^2))*sqrt(sum(b^2)))

option2 <- function(dat, bat, lab) {
  ubat <- levels(bat); N <- ncol(dat)
  alpha <- (rowMeans(dat[,lab==1,drop=FALSE]) + rowMeans(dat[,lab==0,drop=FALSE])) / 2
  n1b <- sapply(ubat, function(b) sum(bat==b & lab==1))
  n0b <- sapply(ubat, function(b) sum(bat==b & lab==0))
  vl <- sapply(seq_along(ubat), function(k) {
    b <- ubat[k]; i <- bat==b
    v1 <- if(n1b[k]>1) rowVars(dat[,i & lab==1,drop=FALSE]) else rep(1e-6, nrow(dat))
    v0 <- if(n0b[k]>1) rowVars(dat[,i & lab==0,drop=FALSE]) else rep(1e-6, nrow(dat))
    pmax((pmax(n1b[k]-1,1)*v1 + pmax(n0b[k]-1,1)*v0) / pmax(n1b[k]+n0b[k]-2,1), 1e-8)
  })
  vp <- rowMeans(vl)
  out <- dat
  for(k in seq_along(ubat)) {
    b <- ubat[k]; i <- bat==b
    m1 <- rowMeans(dat[,i & lab==1,drop=FALSE]); m0 <- rowMeans(dat[,i & lab==0,drop=FALSE])
    gamma_k <- (m1+m0)/2 - alpha
    sf <- sqrt(vp / pmax(vl[,k], 1e-8))
    out[,i] <- sweep(dat[,i,drop=FALSE] - alpha - gamma_k, 1, sf, "*") + alpha
  }
  list(ref=out, vl=vl, vp=vp, ubat=ubat)
}

cohens_d <- function(ref, bat, lab, ubat) {
  sapply(ubat, function(b) {
    i <- bat==b
    n1 <- sum(i & lab==1); n0 <- sum(i & lab==0)
    if(n1<2 || n0<2) return(rep(NA, nrow(ref)))
    m1 <- rowMeans(ref[,i & lab==1,drop=FALSE]); m0 <- rowMeans(ref[,i & lab==0,drop=FALSE])
    v1 <- rowVars(ref[,i & lab==1,drop=FALSE]); v0 <- rowVars(ref[,i & lab==0,drop=FALSE])
    vpool <- ((n1-1)*v1 + (n0-1)*v0) / (n1+n0-2)
    (m1-m0) / sqrt(pmax(vpool, 1e-10))
  })
}  # returns G x B matrix

h3_rows <- list(); h4_rows <- list()

for(nn in 2:5) for(tst in ALL_STUDIES) {
  refs <- train_studies(nn, tst)
  cg <- Reduce(intersect, lapply(c(refs,tst), function(s) rownames(dat_lst[[s]])))
  dat <- log_safe(do.call(cbind, lapply(refs, function(s) dat_lst[[s]][cg,])))
  bat <- as.factor(do.call(c, lapply(refs, function(s) rep(s, ncol(dat_lst[[s]])))))
  lab <- do.call(c, lapply(refs, function(s) bin(label_lst[[s]])))
  ubat <- levels(bat)

  o2 <- option2(dat, bat, lab)
  ref_mo <- suppressMessages(ComBat(dat, bat, mod=model.matrix(~lab), mean.only=TRUE))
  ref_sup<- suppressMessages(ComBat(dat, bat, mod=model.matrix(~lab)))

  # ── H3: Cohen's d preservation ──────────────────────────────────────────────
  d_raw  <- cohens_d(dat,       bat, lab, ubat)
  d_opt2 <- cohens_d(o2$ref,   bat, lab, ubat)
  d_mo   <- cohens_d(ref_mo,   bat, lab, ubat)
  d_sup  <- cohens_d(ref_sup,  bat, lab, ubat)

  # Ratio d_after / d_before per gene per batch (should be ~1 for H3)
  ratio_opt2 <- d_opt2 / (d_raw + 1e-8); ratio_opt2[!is.finite(ratio_opt2)] <- NA
  ratio_mo   <- d_mo   / (d_raw + 1e-8); ratio_mo[!is.finite(ratio_mo)] <- NA
  ratio_sup  <- d_sup  / (d_raw + 1e-8); ratio_sup[!is.finite(ratio_sup)] <- NA

  h3_rows[[length(h3_rows)+1]] <- data.frame(
    n=nn, test=tst,
    # For H3: median ratio, sd ratio (perfect preservation: median~1, sd~0)
    opt2_median_ratio = median(ratio_opt2, na.rm=TRUE),
    opt2_sd_ratio     = sd(ratio_opt2, na.rm=TRUE),
    mo_median_ratio   = median(ratio_mo, na.rm=TRUE),
    mo_sd_ratio       = sd(ratio_mo, na.rm=TRUE),
    sup_median_ratio  = median(ratio_sup, na.rm=TRUE),
    sup_sd_ratio      = sd(ratio_sup, na.rm=TRUE),
    # Absolute deviation: mean(|d_after - d_before|)
    opt2_mad  = mean(abs(d_opt2 - d_raw), na.rm=TRUE),
    mo_mad    = mean(abs(d_mo   - d_raw), na.rm=TRUE),
    sup_mad   = mean(abs(d_sup  - d_raw), na.rm=TRUE),
    stringsAsFactors=FALSE)

  # ── H4: Cosine similarity of pooled class gap vectors ───────────────────────
  gap_raw  <- rowMeans(dat[,lab==1,drop=FALSE])    - rowMeans(dat[,lab==0,drop=FALSE])
  gap_opt2 <- rowMeans(o2$ref[,lab==1,drop=FALSE]) - rowMeans(o2$ref[,lab==0,drop=FALSE])
  gap_mo   <- rowMeans(ref_mo[,lab==1,drop=FALSE])  - rowMeans(ref_mo[,lab==0,drop=FALSE])
  gap_sup  <- rowMeans(ref_sup[,lab==1,drop=FALSE]) - rowMeans(ref_sup[,lab==0,drop=FALSE])

  h4_rows[[length(h4_rows)+1]] <- data.frame(
    n=nn, test=tst,
    cos_opt2_mo   = cosim(gap_opt2, gap_mo),
    cos_opt2_raw  = cosim(gap_opt2, gap_raw),
    cos_mo_raw    = cosim(gap_mo,   gap_raw),
    cos_opt2_sup  = cosim(gap_opt2, gap_sup),
    cos_mo_sup    = cosim(gap_mo,   gap_sup),
    # Magnitude ratio: opt2/mo, mo/raw
    mag_ratio_opt2_mo = sqrt(sum(gap_opt2^2)) / sqrt(sum(gap_mo^2)),
    stringsAsFactors=FALSE)

  cat("done n=",nn," test=",tst,"\n",sep="")
}

h3 <- do.call(rbind, h3_rows)
h4 <- do.call(rbind, h4_rows)
write.csv(h3, file.path(OUT_DIR,"h3_cohens_d.csv"), row.names=FALSE)
write.csv(h4, file.path(OUT_DIR,"h4_gap_similarity.csv"), row.names=FALSE)

cat("\n=== H3: Cohen's d preservation (d_after / d_before, want ~1 for option2) ===\n")
cat("Method   median_ratio  sd_ratio  mean_abs_deviation\n")
cat("opt2    ", round(median(h3$opt2_median_ratio,na.rm=TRUE),3),
    "       ", round(median(h3$opt2_sd_ratio,na.rm=TRUE),3),
    "       ", round(mean(h3$opt2_mad,na.rm=TRUE),3), "\n")
cat("mean.only", round(median(h3$mo_median_ratio,na.rm=TRUE),3),
    "       ", round(median(h3$mo_sd_ratio,na.rm=TRUE),3),
    "       ", round(mean(h3$mo_mad,na.rm=TRUE),3), "\n")
cat("sup     ", round(median(h3$sup_median_ratio,na.rm=TRUE),3),
    "       ", round(median(h3$sup_sd_ratio,na.rm=TRUE),3),
    "       ", round(mean(h3$sup_mad,na.rm=TRUE),3), "\n")

cat("\n=== H4: Cosine similarity of pooled class gap vectors ===\n")
cat("Pair              median cos\n")
cat("opt2 vs mo       ", round(median(h4$cos_opt2_mo),3), "\n")
cat("opt2 vs raw      ", round(median(h4$cos_opt2_raw),3), "\n")
cat("mo vs raw        ", round(median(h4$cos_mo_raw),3), "\n")
cat("opt2 vs sup      ", round(median(h4$cos_opt2_sup),3), "\n")
cat("mo vs sup        ", round(median(h4$cos_mo_sup),3), "\n")
cat("-> h3_cohens_d.csv, h4_gap_similarity.csv\n")
