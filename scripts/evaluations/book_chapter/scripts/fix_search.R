#!/usr/bin/env Rscript
# fix_search.R — find the gene-level cause and a minimal ComBat change that fixes
# the supervised x KNN failure. We vary ONLY step 1, keep step 2 label-blind, and
# measure both a pure-distance metric (nearest-centroid) and the production KNN.
#
# step-1 variants:
#   unsup        ComBat(mod=NULL)                      (baseline, works)
#   sup          ComBat(mod=~class)                    (standard supervised, breaks KNN)
#   sup_meanonly ComBat(mod=~class, mean.only=TRUE)    (kills per-gene batch SCALE delta)
#   sup_nat      ComBat_nat(mod=~class)                (reconstruct w/ batch-only var.pooled)
# Hypothesis: the per-gene batch scale correction (delta.star), interacting with the
# class re-injection, flips the class gap on low-within-variance genes -> mean.only fixes.

suppressMessages({ library(sva); library(matrixStats); library(class) })
source("scripts/ComBat_nat.R")
load("data/TB_real_data.RData")
ALL_STUDIES <- c("GSE37250_SA","USA","India","GSE37250_M","Africa","GSE39941_M")
bin <- function(l) as.integer(ifelse(l %in% c("1",1,"Active"),1L,0L))
train_studies <- function(n,test) ALL_STUDIES[ALL_STUDIES!=test][seq_len(n)]
log_safe <- function(m){ if(max(m,na.rm=TRUE)>100){mn<-min(m,na.rm=TRUE);if(mn<0)m<-m-mn;m<-log2(m+1)}; m }
wvec   <- function(m,l) rowMeans(m[,l==1,drop=FALSE]) - rowMeans(m[,l==0,drop=FALSE])
cosine <- function(a,b) sum(a*b)/sqrt(sum(a^2)*sum(b^2))

step1 <- function(dat,bat,lab,variant){
  if(length(unique(bat))<2) return(dat)
  switch(variant,
    unsup        = suppressMessages(ComBat(dat,batch=bat,mod=NULL)),
    sup          = suppressMessages(ComBat(dat,batch=bat,mod=model.matrix(~lab))),
    sup_meanonly = suppressMessages(ComBat(dat,batch=bat,mod=model.matrix(~lab),mean.only=TRUE)),
    sup_nat      = suppressMessages(ComBat_nat(dat,batch=bat,mod=model.matrix(~lab))))
}
align_test <- function(ref,tst){           # label-blind step 2 (same for all variants)
  comb<-cbind(ref,tst); cb<-c(rep(1L,ncol(ref)),rep(2L,ncol(tst)))
  out<-suppressMessages(ComBat(comb,batch=cb,mod=NULL,ref.batch=1L))
  out[,(ncol(ref)+1):ncol(out),drop=FALSE]
}
cen_bacc <- function(ref,tgt,rl,tl){
  ca<-rowMeans(ref[,rl==1,drop=FALSE]); ch<-rowMeans(ref[,rl==0,drop=FALSE])
  da<-sqrt(colSums((tgt-ca)^2)); dh<-sqrt(colSums((tgt-ch)^2)); pred<-as.integer(da<dh)
  tp<-sum(pred==1&tl==1);tn<-sum(pred==0&tl==0);fp<-sum(pred==1&tl==0);fn<-sum(pred==0&tl==1)
  mean(c(if((tp+fn)>0)tp/(tp+fn)else NA,if((tn+fp)>0)tn/(tn+fp)else NA),na.rm=TRUE)
}
# faithful production KNN (k chosen by the same CV as predKNN_pp), test prediction
knn_prod <- function(ref,tgt,rl,tl,n_genes=1000){
  cg<-intersect(rownames(ref),rownames(tgt)); ref<-ref[cg,];tgt<-tgt[cg,]
  hvg<-order(rowVars(ref),decreasing=TRUE)[seq_len(min(n_genes,nrow(ref)))]
  X<-t(ref[hvg,,drop=FALSE]); Y<-t(tgt[hvg,,drop=FALSE]); ns<-nrow(X)
  ks<-c(3,5,7,9,11); ks<-ks[ks<ns]; if(!length(ks)) ks<-min(3,ns-1)
  acc<-sapply(ks,function(k){ f<-cut(seq_len(ns),5,labels=FALSE); ok<-0;tot<-0
    for(fo in 1:5){ te<-which(f==fo);tr<-which(f!=fo); if(length(te)&&length(tr)){
      p<-class::knn(X[tr,,drop=FALSE],X[te,,drop=FALSE],rl[tr],k=k); ok<-ok+sum(p==rl[te]);tot<-tot+length(te)}}
    if(tot)ok/tot else 0})
  k<-ks[which.max(acc)]
  pred<-as.integer(as.character(class::knn(X,Y,rl,k=k)))
  tp<-sum(pred==1&tl==1);tn<-sum(pred==0&tl==0);fp<-sum(pred==1&tl==0);fn<-sum(pred==0&tl==1)
  den<-sqrt(as.double(tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
  list(bacc=mean(c(if((tp+fn)>0)tp/(tp+fn)else NA,if((tn+fp)>0)tn/(tn+fp)else NA),na.rm=TRUE),
       mcc=if(den>0)(tp*tn-fp*fn)/den else 0)
}

variants<-c("unsup","sup","sup_meanonly","sup_nat")
rows<-list()
for(n in 2:5) for(test in ALL_STUDIES){
  refs<-train_studies(n,test); cg<-Reduce(intersect,lapply(c(refs,test),function(s)rownames(dat_lst[[s]])))
  dat<-log_safe(do.call(cbind,lapply(refs,function(s)dat_lst[[s]][cg,])))
  bat<-do.call(c,lapply(refs,function(s)rep(s,ncol(dat_lst[[s]])))); lab<-do.call(c,lapply(refs,function(s)bin(label_lst[[s]])))
  tst<-log_safe(dat_lst[[test]][cg,]); ltst<-bin(label_lst[[test]]); w_test<-wvec(tst,ltst)
  for(v in variants){
    ref<-tryCatch(step1(dat,bat,lab,v),error=function(e)NULL); if(is.null(ref))next
    tgt<-align_test(ref,tst); kp<-knn_prod(ref,tgt,lab,ltst)
    rows[[length(rows)+1]]<-data.frame(n=n,test=test,variant=v,
      cos_test=cosine(wvec(ref,lab),w_test), cen_bacc=cen_bacc(ref,tgt,lab,ltst),
      knn_bacc=kp$bacc, knn_mcc=kp$mcc)
  }
  cat(sprintf("  done n=%d test=%s\n",n,test))
}
R<-do.call(rbind,rows); write.csv(R,"outputs/diagnostics/hypothesis_tests/fix_search.csv",row.names=FALSE)
cat("\n── Medians by step-1 variant (24 scenarios) ──\n")
agg<-do.call(rbind,lapply(variants,function(v){s<-R[R$variant==v,]
  data.frame(variant=v, cos_test=round(median(s$cos_test),3), cen_bacc=round(median(s$cen_bacc),3),
             knn_bacc=round(median(s$knn_bacc),3), knn_mcc=round(median(s$knn_mcc),3),
             n_knn_below_chance=sum(s$knn_mcc<0))}))
print(agg,row.names=FALSE)
cat("\n-> fix_search.csv\n")
