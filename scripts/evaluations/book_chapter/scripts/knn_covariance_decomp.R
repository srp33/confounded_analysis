#!/usr/bin/env Rscript
# knn_covariance_decomp.R
# Decompose the KNN (raw-Euclidean) class signal by VARIANCE direction, per the
# covariance angle: KNN distance d^2 = sum_k (proj on PC_k)^2, and high-variance
# PCs dominate which points are neighbours. Question: under supervised ComBat, does
# the class signal that KNN actually weights (high-variance PCs) point the WRONG way
# relative to the held-out test, while the right signal hides in low-variance PCs
# (so whitening / LDA recovers it)?
#
# For PC_k of the step-1 reference:  lambda_k = variance,
#   a_k = training class-mean gap on PC_k,   t_k = TEST class-mean gap on PC_k.
# raw (KNN) alignment       = sum_k a_k t_k          (basis-free = <w_train,w_test>)
# whitened (LDA) alignment  = sum_k a_k t_k / lambda_k
# We show WHICH variance bands carry each, for a FAILING and a WORKING scenario.

suppressMessages({ library(sva); library(matrixStats); library(ggplot2); library(dplyr) })
load("data/TB_real_data.RData"); OUT_DIR<-"outputs/diagnostics/hypothesis_tests"
ALL_STUDIES<-c("GSE37250_SA","USA","India","GSE37250_M","Africa","GSE39941_M")
bin<-function(l) as.integer(ifelse(l%in%c("1",1,"Active"),1L,0L))
train_studies<-function(n,test) ALL_STUDIES[ALL_STUDIES!=test][seq_len(n)]
log_safe<-function(m){if(max(m,na.rm=TRUE)>100){mn<-min(m,na.rm=TRUE);if(mn<0)m<-m-mn;m<-log2(m+1)};m}
step1<-function(dat,bat,lab,sup) if(length(unique(bat))>=2)
  suppressMessages(ComBat(dat,batch=bat,mod=if(sup)model.matrix(~lab)else NULL)) else dat

decomp<-function(n,test,sup,K=80){
  refs<-train_studies(n,test); cg<-Reduce(intersect,lapply(c(refs,test),function(s)rownames(dat_lst[[s]])))
  dat<-log_safe(do.call(cbind,lapply(refs,function(s)dat_lst[[s]][cg,])))
  bat<-do.call(c,lapply(refs,function(s)rep(s,ncol(dat_lst[[s]]))))
  lab<-do.call(c,lapply(refs,function(s)bin(label_lst[[s]])))
  tst<-log_safe(dat_lst[[test]][cg,]); ltst<-bin(label_lst[[test]])
  ref<-step1(dat,bat,lab,sup)
  K<-min(K,ncol(ref)-2); pca<-prcomp(t(ref),rank.=K)
  Rp<-pca$x; Tp<-predict(pca,t(tst))
  lambda<-pca$sdev[1:K]^2
  a<-colMeans(Rp[lab==1,,drop=FALSE])-colMeans(Rp[lab==0,,drop=FALSE])   # train class gap per PC
  t_<-colMeans(Tp[ltst==1,,drop=FALSE])-colMeans(Tp[ltst==0,,drop=FALSE])# test  class gap per PC
  data.frame(scen=sprintf("n%d_%s",n,test), sup=sup, pc=1:K, lambda=lambda,
             a=a, t=t_, contrib_raw=a*t_, contrib_white=a*t_/lambda)
}

scn<-list(list(5,"GSE37250_SA","FAILING"), list(3,"USA","WORKING"))
allrows<-list()
for(s in scn) for(sup in c(FALSE,TRUE)){
  d<-decomp(s[[1]],s[[2]],sup); d$label<-s[[3]]; allrows[[length(allrows)+1]]<-d
}
D<-do.call(rbind,allrows)
write.csv(D,file.path(OUT_DIR,"knn_covariance_decomp.csv"),row.names=FALSE)

# ── summary: where does the raw (KNN) alignment come from? ────────────────────
cat("Raw (KNN) alignment = sum a_k*t_k ; whitened (LDA) = sum a_k*t_k/lambda_k\n")
cat("Fraction of |raw| contributed by the top-5 highest-variance PCs, and its SIGN:\n\n")
for(s in scn) for(sup in c(FALSE,TRUE)){
  x<-D[D$scen==sprintf("n%d_%s",s[[1]],s[[2]]) & D$sup==sup,]
  x<-x[order(-x$lambda),]
  raw_tot<-sum(x$contrib_raw); top5<-sum(x$contrib_raw[1:5]); wh_tot<-sum(x$contrib_white)
  cat(sprintf("  %-14s %-12s raw=%+8.2f (top5 PCs=%+8.2f, %.0f%% of |raw|, sign %s) | whitened=%+8.3f\n",
      s[[3]], ifelse(sup,"supervised","unsupervised"), raw_tot, top5,
      100*abs(top5)/sum(abs(x$contrib_raw)), ifelse(top5<0,"WRONG","right"), wh_tot))
}

# ── figure: per-PC signed contribution to KNN alignment vs PC variance ────────
D$step1<-factor(ifelse(D$sup,"Supervised","Unsupervised"),levels=c("Unsupervised","Supervised"))
D$label<-factor(D$label,levels=c("FAILING","WORKING"))
# cumulative raw-alignment as KNN accumulates PCs from HIGH to LOW variance, normalised
Dc <- D %>% group_by(label, step1) %>% arrange(desc(lambda), .by_group=TRUE) %>%
  mutate(rank = row_number(), cum = cumsum(contrib_raw), cum_norm = cum / sum(abs(contrib_raw)))
p<-ggplot(Dc, aes(rank, cum_norm, colour=step1)) +
  geom_hline(yintercept=0,linewidth=0.3,colour="grey30")+
  geom_line(linewidth=1)+
  facet_wrap(~label)+
  scale_colour_manual(values=c(Unsupervised="#4e79a7",Supervised="#e15759"),name=NULL)+
  labs(title="Where KNN's class signal comes from, by variance direction",
       subtitle="Cumulative train-test class alignment as PCs are added HIGHEST-variance first.\nKNN is dominated by the left side. Below zero = neighbours of the WRONG class.",
       x="Number of PCs included (ordered high -> low variance)",
       y="Cumulative class alignment (normalised)")+
  theme_minimal(base_size=12)+ theme(legend.position=c(0.5,0.5))
ggsave(file.path(OUT_DIR,"knn_covariance_decomp.png"),p,width=11,height=5.5,dpi=200)
ggsave(file.path(OUT_DIR,"knn_covariance_decomp.pdf"),p,width=11,height=5.5)

# ── OOD check: is the TEST out-of-distribution along high-variance dirs? ───────
cat("\nOOD check — median |test PC score| / sqrt(lambda) (1 ~ in-distribution):\n")
for(s in scn) for(sup in c(FALSE,TRUE)){
  refs<-train_studies(s[[1]],s[[2]]); cg<-Reduce(intersect,lapply(c(refs,s[[2]]),function(z)rownames(dat_lst[[z]])))
  dat<-log_safe(do.call(cbind,lapply(refs,function(z)dat_lst[[z]][cg,])))
  bat<-do.call(c,lapply(refs,function(z)rep(z,ncol(dat_lst[[z]])))); lab<-do.call(c,lapply(refs,function(z)bin(label_lst[[z]])))
  tst<-log_safe(dat_lst[[s[[2]]]][cg,]); ref<-step1(dat,bat,lab,sup)
  K<-min(20,ncol(ref)-2); pca<-prcomp(t(ref),rank.=K); Tp<-predict(pca,t(tst))
  z<-sweep(abs(Tp),2,pca$sdev[1:K],"/")
  cat(sprintf("  %-14s %-12s top5 PCs=%.2f  PCs6-20=%.2f\n", s[[3]], ifelse(sup,"supervised","unsupervised"),
      median(z[,1:5]), median(z[,6:K])))
}
cat("\n-> knn_covariance_decomp.csv/.png\n")
