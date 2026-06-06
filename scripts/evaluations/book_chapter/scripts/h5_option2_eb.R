#!/usr/bin/env Rscript
# h5_option2_eb.R — Does EB-shrunk within-class variance improve Option2?
# H5: batches with small n have noisy vl estimates; EB shrinkage stabilizes them.
# Tests option2_eb vs option2 vs mean.only across 24 scenarios.

suppressMessages({ library(sva); library(matrixStats); library(class) })
setwd("/home/phr23/confounded_analysis/scripts/evaluations/book_chapter")
load("data/TB_real_data.RData")
OUT_DIR <- "outputs/diagnostics/hypothesis_tests"

ALL_STUDIES <- c("GSE37250_SA","USA","India","GSE37250_M","Africa","GSE39941_M")
bin <- function(l) as.integer(ifelse(l %in% c("1",1,"Active"), 1L, 0L))
train_studies <- function(n, test) ALL_STUDIES[ALL_STUDIES!=test][seq_len(n)]
log_safe <- function(m) { if(max(m,na.rm=TRUE)>100){ mn<-min(m,na.rm=TRUE); if(mn<0) m<-m-mn; m<-log2(m+1) }; m }

aprior <- function(gg) { m<-mean(gg); s2<-var(gg); (2*s2+m^2)/s2 }
bprior <- function(gg) { m<-mean(gg); s2<-var(gg); (m*s2+m^3)/s2 }
postvar_fn <- function(sum2, n, a, b) (0.5*sum2+b)/(n/2+a-1)

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
  out
}

option2_eb <- function(dat, bat, lab) {
  ubat <- levels(bat); N <- ncol(dat)
  alpha <- (rowMeans(dat[,lab==1,drop=FALSE]) + rowMeans(dat[,lab==0,drop=FALSE])) / 2
  n1b <- sapply(ubat, function(b) sum(bat==b & lab==1))
  n0b <- sapply(ubat, function(b) sum(bat==b & lab==0))
  vl_raw <- sapply(seq_along(ubat), function(k) {
    b <- ubat[k]; i <- bat==b
    v1 <- if(n1b[k]>1) rowVars(dat[,i & lab==1,drop=FALSE]) else rep(1e-6, nrow(dat))
    v0 <- if(n0b[k]>1) rowVars(dat[,i & lab==0,drop=FALSE]) else rep(1e-6, nrow(dat))
    pmax((pmax(n1b[k]-1,1)*v1 + pmax(n0b[k]-1,1)*v0) / pmax(n1b[k]+n0b[k]-2,1), 1e-8)
  })
  vp <- rowMeans(vl_raw)
  # EB shrinkage of within-class variances using inverse-gamma prior
  # Prior fitted per-batch from gene-level distribution of vl_raw
  vl_star <- vl_raw
  for(k in seq_along(ubat)) {
    df_b <- n1b[k] + n0b[k] - 2
    if(df_b < 1) next
    vl_b <- vl_raw[,k]
    if(var(vl_b) > 1e-12) {
      a_b <- aprior(vl_b); s_b <- bprior(vl_b)
      # Posterior mode: treats vl_b as (sum_sq / df_b), so sum_sq = vl_b * df_b
      vs <- postvar_fn(vl_b * df_b, df_b, a_b, s_b)
      vl_star[,k] <- pmax(vs, 1e-8)
    }
  }
  out <- dat
  for(k in seq_along(ubat)) {
    b <- ubat[k]; i <- bat==b
    m1 <- rowMeans(dat[,i & lab==1,drop=FALSE]); m0 <- rowMeans(dat[,i & lab==0,drop=FALSE])
    gamma_k <- (m1+m0)/2 - alpha
    sf <- sqrt(vp / pmax(vl_star[,k], 1e-8))
    out[,i] <- sweep(dat[,i,drop=FALSE] - alpha - gamma_k, 1, sf, "*") + alpha
  }
  out
}

align_test <- function(ref, test) {
  comb <- cbind(ref, test); cb <- c(rep(1L, ncol(ref)), rep(2L, ncol(test)))
  suppressMessages(ComBat(comb, batch=cb, mod=NULL, ref.batch=1L))[,(ncol(ref)+1):ncol(comb),drop=FALSE]
}

knn_mcc <- function(train, test_cor, lab, lab_test, k=5) {
  hvg <- order(rowVars(train), decreasing=TRUE)[1:min(1000,nrow(train))]
  pred <- as.integer(as.character(knn(t(train[hvg,]), t(test_cor[hvg,]), lab, k=k)))
  tp <- sum(pred==1 & lab_test==1); fp <- sum(pred==1 & lab_test==0)
  tn <- sum(pred==0 & lab_test==0); fn <- sum(pred==0 & lab_test==1)
  denom <- sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
  if(denom==0) 0 else (tp*tn - fp*fn) / denom
}

rows <- list()
for(nn in 2:5) for(tst in ALL_STUDIES) {
  refs <- train_studies(nn, tst)
  cg <- Reduce(intersect, lapply(c(refs,tst), function(s) rownames(dat_lst[[s]])))
  dat  <- log_safe(do.call(cbind, lapply(refs, function(s) dat_lst[[s]][cg,])))
  tdat <- log_safe(dat_lst[[tst]][cg,])
  bat  <- as.factor(do.call(c, lapply(refs, function(s) rep(s, ncol(dat_lst[[s]])))))
  lab  <- do.call(c, lapply(refs, function(s) bin(label_lst[[s]])))
  ltst <- bin(label_lst[[tst]])

  ref_combat  <- suppressMessages(ComBat(dat, bat, mod=NULL))
  ref_mo      <- suppressMessages(ComBat(dat, bat, mod=model.matrix(~lab), mean.only=TRUE))
  ref_o2      <- option2(dat, bat, lab)
  ref_o2eb    <- option2_eb(dat, bat, lab)

  # Align test set to each reference (unsupervised step 2)
  tst_combat  <- align_test(ref_combat,  tdat)
  tst_mo      <- align_test(ref_mo,      tdat)
  tst_o2      <- align_test(ref_o2,      tdat)
  tst_o2eb    <- align_test(ref_o2eb,    tdat)

  rows[[length(rows)+1]] <- data.frame(
    n=nn, test=tst,
    mcc_combat  = knn_mcc(ref_combat,  tst_combat,  lab, ltst),
    mcc_mo      = knn_mcc(ref_mo,      tst_mo,      lab, ltst),
    mcc_opt2    = knn_mcc(ref_o2,      tst_o2,      lab, ltst),
    mcc_opt2_eb = knn_mcc(ref_o2eb,    tst_o2eb,    lab, ltst),
    stringsAsFactors=FALSE)
  cat("done n=",nn," test=",tst,"\n",sep="")
}
res <- do.call(rbind, rows)
write.csv(res, file.path(OUT_DIR,"h5_option2_eb_results.csv"), row.names=FALSE)

cat("\n=== H5: KNN MCC across 24 scenarios ===\n")
cat(sprintf("%-12s  median  n_below_chance\n","method"))
for(col in c("mcc_combat","mcc_mo","mcc_opt2","mcc_opt2_eb")) {
  x <- res[[col]]
  cat(sprintf("%-12s  %6.3f  %d\n", col, median(x), sum(x<0)))
}
cat("-> h5_option2_eb_results.csv\n")
