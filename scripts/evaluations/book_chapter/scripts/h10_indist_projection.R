#!/usr/bin/env Rscript
# h10_indist_projection.R
# Question: is held-out OOD required for the KNN collapse, or is label-blind
# STEP-2 PROJECTION onto the distorted reference sufficient — even for
# IN-DISTRIBUTION samples (drawn from a study that helped build the reference)?
#
# Three coordinates for the SAME geometry, contrasted:
#   (1) resub        : adjusted training samples, placed by the label-aware joint fit (knn.cv LOO)
#   (2) indist_proj  : samples from a TRAINING study, held out of the fit, then
#                      re-projected onto the reference by unsupervised step-2
#   (3) ood_proj     : a held-out STUDY, projected onto the reference by step-2
# If (2) collapses like (3) while (1) survives -> OOD is NOT required; the
# label-blind projection onto a distorted reference is the operative cause.

suppressMessages({ library(sva); library(matrixStats); library(class) })
setwd("/home/phr23/confounded_analysis/scripts/evaluations/book_chapter")
load("data/TB_real_data.RData")
OUT_DIR <- "outputs/diagnostics/hypothesis_tests"

ALL_STUDIES <- c("GSE37250_SA","USA","India","GSE37250_M","Africa","GSE39941_M")
bin <- function(l) as.integer(ifelse(l %in% c("1",1,"Active"), 1L, 0L))
train_studies <- function(n, test) ALL_STUDIES[ALL_STUDIES!=test][seq_len(n)]
log_safe <- function(m) { if(max(m,na.rm=TRUE)>100){mn<-min(m,na.rm=TRUE);if(mn<0)m<-m-mn;m<-log2(m+1)};m }

mcc <- function(pred, truth) {
  tp <- as.numeric(sum(pred==1 & truth==1)); fp <- as.numeric(sum(pred==1 & truth==0))
  tn <- as.numeric(sum(pred==0 & truth==0)); fn <- as.numeric(sum(pred==0 & truth==1))
  d <- sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
  if(is.na(d) || d==0) 0 else (tp*tn-fp*fn)/d
}
align_test <- function(ref, test) {
  comb <- cbind(ref, test); cb <- c(rep(1L,ncol(ref)), rep(2L,ncol(test)))
  suppressMessages(ComBat(comb, batch=cb, mod=NULL, ref.batch=1L))[,(ncol(ref)+1):ncol(comb),drop=FALSE]
}
hvg_knn <- function(train, query, ltrain, k=5) {
  hv <- order(rowVars(train), decreasing=TRUE)[1:min(1000,nrow(train))]
  as.integer(as.character(knn(t(train[hv,,drop=FALSE]), t(query[hv,,drop=FALSE]), ltrain, k=k)))
}

set.seed(11)
rows <- list()
n <- 5
for(tst in ALL_STUDIES) {
  refs <- train_studies(n, tst)
  cg   <- Reduce(intersect, lapply(c(refs,tst), function(s) rownames(dat_lst[[s]])))
  tdat <- log_safe(dat_lst[[tst]][cg,]); ltst <- bin(label_lst[[tst]])

  # raw per-study matrices/labels
  raw <- lapply(refs, function(s) log_safe(dat_lst[[s]][cg,]))
  names(raw) <- refs
  lb  <- lapply(refs, function(s) bin(label_lst[[s]])); names(lb) <- refs

  # ----- per training study: SAME withheld samples placed two ways -----
  # For study s, split 60/40. The "build" set = other studies + s_build.
  # The SAME s_held samples are then placed by:
  #   (a) STEP-1 resub : include s_held in the supervised joint fit -> label-aware coords
  #   (b) STEP-2 proj  : exclude s_held from the fit, align via unsupervised step-2 -> label-blind
  # Neighbors in both cases = the label-aware reference build coords. Only the
  # placement mechanism for s_held differs. Plus (c) the OOD study via step-2.
  resub_mccs <- c(); step1_held_mccs <- c(); step2_held_mccs <- c(); ood_mccs <- c()
  for(s in refs) {
    ls <- lb[[s]]; ns <- length(ls)
    if(sum(ls==1) < 4 || sum(ls==0) < 4) next  # need class balance to split
    idx1 <- which(ls==1); idx0 <- which(ls==0)
    b1 <- sample(idx1, ceiling(0.6*length(idx1))); b0 <- sample(idx0, ceiling(0.6*length(idx0)))
    build_idx <- c(b1,b0); held_idx <- setdiff(seq_len(ns), build_idx)
    if(length(held_idx) < 4 || length(unique(ls[held_idx])) < 2) next

    others <- setdiff(refs, s)
    oth_dat <- do.call(cbind, raw[others])
    oth_bat <- do.call(c, lapply(others, function(o) rep(o, ncol(raw[[o]]))))
    oth_lab <- do.call(c, lb[others])
    s_build <- raw[[s]][,build_idx,drop=FALSE]; lab_build <- ls[build_idx]
    s_held  <- raw[[s]][,held_idx,drop=FALSE];  lab_held  <- ls[held_idx]

    # Reference WITHOUT s_held (s still in-distribution via s_build)
    dat_wo <- cbind(oth_dat, s_build)
    bat_wo <- as.factor(c(oth_bat, rep(s, ncol(s_build))))
    lab_wo <- c(oth_lab, lab_build)
    ref_wo <- suppressMessages(ComBat(dat_wo, bat_wo, mod=model.matrix(~lab_wo)))

    # Reference WITH s_held in the step-1 fit (same s batch, just more samples)
    dat_wi <- cbind(oth_dat, raw[[s]])
    bat_wi <- as.factor(c(oth_bat, rep(s, ncol(raw[[s]]))))
    lab_wi <- c(oth_lab, ls)
    ref_wi <- suppressMessages(ComBat(dat_wi, bat_wi, mod=model.matrix(~lab_wi)))
    # columns of ref_wi corresponding to s_held vs the rest
    s_cols <- (ncol(oth_dat)+1):ncol(ref_wi)        # all of study s within ref_wi
    held_cols <- s_cols[held_idx]                    # the withheld subset
    train_cols<- setdiff(seq_len(ncol(ref_wi)), held_cols)

    # (a) STEP-1 resub: s_held placed by joint fit (cols of ref_wi), neighbors = everything else
    pa <- hvg_knn(ref_wi[,train_cols,drop=FALSE], ref_wi[,held_cols,drop=FALSE], factor(lab_wi[train_cols]))
    step1_held_mccs <- c(step1_held_mccs, mcc(pa, lab_held))

    # (b) STEP-2 proj: same s_held placed by unsupervised step-2 onto ref_wo
    held_proj <- align_test(ref_wo, s_held)
    pb <- hvg_knn(ref_wo, held_proj, factor(lab_wo))
    step2_held_mccs <- c(step2_held_mccs, mcc(pb, lab_held))

    # (c) OOD study via step-2 onto ref_wo, for scale
    ood_proj <- align_test(ref_wo, tdat)
    pc <- hvg_knn(ref_wo, ood_proj, factor(lab_wo))
    ood_mccs <- c(ood_mccs, mcc(pc, ltst))

    # plain resub on the with-fit reference (LOO), for the classic baseline
    hv <- order(rowVars(ref_wi), decreasing=TRUE)[1:1000]
    pr <- as.integer(as.character(knn.cv(t(ref_wi[hv,]), factor(lab_wi), k=5)))
    resub_mccs <- c(resub_mccs, mcc(pr, lab_wi))
  }

  rows[[length(rows)+1]] <- data.frame(
    test=tst,
    mcc_resub_loo      = mean(resub_mccs, na.rm=TRUE),      # all training, label-aware, LOO
    mcc_step1_held     = mean(step1_held_mccs, na.rm=TRUE), # s_held placed by step-1 joint fit
    mcc_step2_held     = mean(step2_held_mccs, na.rm=TRUE), # SAME s_held placed by step-2 projection
    mcc_ood_step2      = mean(ood_mccs, na.rm=TRUE),        # OOD study via step-2
    n_study_splits     = length(step2_held_mccs),
    stringsAsFactors=FALSE)
  cat("done test=",tst,
      " resub=",round(mean(resub_mccs,na.rm=TRUE),3),
      " step1_held=",round(mean(step1_held_mccs,na.rm=TRUE),3),
      " step2_held=",round(mean(step2_held_mccs,na.rm=TRUE),3),
      " ood=",round(mean(ood_mccs,na.rm=TRUE),3),"\n",sep="")
}
res <- do.call(rbind, rows)
write.csv(res, file.path(OUT_DIR,"h10_indist_projection.csv"), row.names=FALSE)

cat("\n=== H10: Same withheld samples, placed by STEP-1 vs STEP-2 ===\n")
cat(sprintf("%-18s median   mean\n","coordinate"))
for(col in c("mcc_resub_loo","mcc_step1_held","mcc_step2_held","mcc_ood_step2")) {
  x <- res[[col]]; cat(sprintf("%-18s %6.3f  %6.3f\n", col, median(x,na.rm=TRUE), mean(x,na.rm=TRUE)))
}
cat("\nstep1_held vs step2_held are the SAME in-distribution samples; only placement differs.\n")
cat("If step2_held collapses to ~ood_step2 while step1_held survives,\n")
cat("then label-blind STEP-2 projection is the operative cause -- OOD is NOT required.\n")
cat("-> h10_indist_projection.csv\n")
