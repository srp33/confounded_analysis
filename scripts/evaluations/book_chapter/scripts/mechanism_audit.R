#!/usr/bin/env Rscript
# mechanism_audit.R
# Stress-tests the mechanism claims against the critiques:
#  (1) Can ComBat flip a gene's class-difference sign?  -> count sign flips.
#  (5) Does the per-gene class direction actually change between unsup/sup step1,
#      and relative to the RAW data and the TEST's own class direction? -> cosines.
#  (2/4) KNN cares about DISTANCES, not a projection axis. Measure a direction-free,
#      pure-distance quantity: is each corrected test sample Euclidean-closer to the
#      correct-class training centroid?  (nearest-centroid accuracy in gene space.)

suppressMessages({ library(sva); library(matrixStats) })
DATA_FILE <- "data/TB_real_data.RData"; OUT_DIR <- "outputs/diagnostics/hypothesis_tests"
load(DATA_FILE)
ALL_STUDIES <- c("GSE37250_SA","USA","India","GSE37250_M","Africa","GSE39941_M")
bin <- function(l) as.integer(ifelse(l %in% c("1",1,"Active"),1L,0L))
train_studies <- function(n,test) ALL_STUDIES[ALL_STUDIES!=test][seq_len(n)]
log_safe <- function(m){ if(max(m,na.rm=TRUE)>100){mn<-min(m,na.rm=TRUE);if(mn<0)m<-m-mn;m<-log2(m+1)}; m }

two_step <- function(dat,bat,lab,tst,ltst,sup1,sup2){
  dat_c <- if(length(unique(bat))>=2)
    suppressMessages(ComBat(dat,batch=bat,mod=if(sup1) model.matrix(~lab) else NULL)) else dat
  comb<-cbind(dat_c,tst); cb<-c(rep(1L,ncol(dat_c)),rep(2L,ncol(tst)))
  mod2<-if(sup2) model.matrix(~c(lab,ltst)) else NULL
  out<-suppressMessages(ComBat(comb,batch=cb,mod=mod2,ref.batch=1L))
  list(ref=dat_c, tgt=out[,(ncol(dat_c)+1):ncol(out),drop=FALSE])
}
wvec  <- function(m,l) rowMeans(m[,l==1,drop=FALSE]) - rowMeans(m[,l==0,drop=FALSE])
cosine<- function(a,b) sum(a*b)/sqrt(sum(a^2)*sum(b^2))

# direction-free, pure-distance class assignment (nearest training centroid)
centroid_acc <- function(ref,tgt,rl,tl,n_genes=1000){
  cg<-intersect(rownames(ref),rownames(tgt)); ref<-ref[cg,];tgt<-tgt[cg,]
  hvg<-order(rowVars(ref),decreasing=TRUE)[seq_len(min(n_genes,nrow(ref)))]
  ref<-ref[hvg,];tgt<-tgt[hvg,]
  ca<-rowMeans(ref[,rl==1,drop=FALSE]); ch<-rowMeans(ref[,rl==0,drop=FALSE])
  da<-sqrt(colSums((tgt-ca)^2)); dh<-sqrt(colSums((tgt-ch)^2))
  pred<-as.integer(da<dh)
  tp<-sum(pred==1&tl==1);tn<-sum(pred==0&tl==0);fp<-sum(pred==1&tl==0);fn<-sum(pred==0&tl==1)
  list(bacc=mean(c(if((tp+fn)>0)tp/(tp+fn) else NA, if((tn+fp)>0)tn/(tn+fp) else NA),na.rm=TRUE),
       frac_act_correct = if(sum(tl==1)>0) mean(da[tl==1]<dh[tl==1]) else NA,
       frac_hlt_correct = if(sum(tl==0)>0) mean(dh[tl==0]<da[tl==0]) else NA)
}

rows<-list()
for(n in 2:5) for(test in ALL_STUDIES){
  refs<-train_studies(n,test)
  cg<-Reduce(intersect,lapply(c(refs,test),function(s) rownames(dat_lst[[s]])))
  dat<-log_safe(do.call(cbind,lapply(refs,function(s) dat_lst[[s]][cg,])))
  bat<-do.call(c,lapply(refs,function(s) rep(s,ncol(dat_lst[[s]]))))
  lab<-do.call(c,lapply(refs,function(s) bin(label_lst[[s]])))
  tst<-log_safe(dat_lst[[test]][cg,]); ltst<-bin(label_lst[[test]])

  ref_u<-two_step(dat,bat,lab,tst,ltst,FALSE,FALSE)   # combat
  ref_s<-two_step(dat,bat,lab,tst,ltst,TRUE ,FALSE)   # sup step1, label-blind step2 (isolates step1)

  # ---- direction metrics (step1 references), all on common genes ----
  g<-Reduce(intersect,list(rownames(dat),rownames(ref_u$ref),rownames(ref_s$ref),rownames(tst)))
  w_raw <- wvec(dat[g,],lab); w_u<-wvec(ref_u$ref[g,],lab); w_s<-wvec(ref_s$ref[g,],lab)
  w_test<- wvec(tst[g,],ltst)
  rows[[length(rows)+1]]<-data.frame(
    n=n,test=test,
    flip_unsup = mean(sign(w_raw)!=sign(w_u)),           # critique 1
    flip_sup   = mean(sign(w_raw)!=sign(w_s)),
    cos_raw_unsup = cosine(w_raw,w_u),                    # critique 5
    cos_raw_sup   = cosine(w_raw,w_s),
    cos_unsup_sup = cosine(w_u,w_s),
    cos_unsup_test= cosine(w_u,w_test),                   # does train axis match TEST's own?
    cos_sup_test  = cosine(w_s,w_test),
    # ---- pure-distance (direction-free) KNN-style metric, critiques 2/4 ----
    cen_bacc_unsup = centroid_acc(ref_u$ref,ref_u$tgt,lab,ltst)$bacc,
    cen_bacc_sup   = centroid_acc(ref_s$ref,ref_s$tgt,lab,ltst)$bacc)
  cat(sprintf("  done n=%d test=%s\n",n,test))
}
A<-do.call(rbind,rows)
write.csv(A,file.path(OUT_DIR,"mechanism_audit.csv"),row.names=FALSE)

f <- function(x) round(median(x,na.rm=TRUE),3)
cat("\n── (1) per-gene class-difference SIGN FLIPS (median frac of genes) ──\n")
cat(sprintf("   raw -> unsup step1: %.4f      raw -> sup step1: %.4f\n", f(A$flip_unsup), f(A$flip_sup)))
cat("\n── (5) direction change (median cosine; 1 = identical direction) ──\n")
cat(sprintf("   cos(raw, unsup)=%.3f   cos(raw, sup)=%.3f   cos(unsup, sup)=%.3f\n",
            f(A$cos_raw_unsup), f(A$cos_raw_sup), f(A$cos_unsup_sup)))
cat(sprintf("   cos(train axis, TEST's own axis): unsup=%.3f   sup=%.3f\n",
            f(A$cos_unsup_test), f(A$cos_sup_test)))
cat("\n── (2/4) DIRECTION-FREE distance test: nearest-centroid balanced acc ──\n")
cat(sprintf("   unsup step1 = %.3f      sup step1 = %.3f\n", f(A$cen_bacc_unsup), f(A$cen_bacc_sup)))
cat("\n  -> mechanism_audit.csv\n")
