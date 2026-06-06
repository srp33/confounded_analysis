#!/usr/bin/env Rscript
# assumption_tests.R  v2 — run on a FAILING scenario and a WORKING one for contrast.
#  #12 Is w_sup the PARTIAL (batch-adjusted) class coefficient, and does the partial
#      differ from the MARGINAL class difference under confounding?
#  #14 Is batch-class confounding CAUSAL?  De-confound (balance class within each
#      study), re-run, see if the failure disappears.
#  #15 Does S_W^-1 whitening (LDA) recover what raw Euclidean (KNN) loses?
# Methodology kept consistent with mechanism_audit.R: two_step(sup1=TRUE, sup2=FALSE)
# isolates step 1 while keeping the test step-2 aligned (label-blind).

suppressMessages({ library(sva); library(matrixStats) })
load("data/TB_real_data.RData")
ALL_STUDIES <- c("GSE37250_SA","USA","India","GSE37250_M","Africa","GSE39941_M")
bin <- function(l) as.integer(ifelse(l %in% c("1",1,"Active"),1L,0L))
train_studies <- function(n,test) ALL_STUDIES[ALL_STUDIES!=test][seq_len(n)]
log_safe <- function(m){ if(max(m,na.rm=TRUE)>100){mn<-min(m,na.rm=TRUE);if(mn<0)m<-m-mn;m<-log2(m+1)}; m }
wvec   <- function(m,l) rowMeans(m[,l==1,drop=FALSE]) - rowMeans(m[,l==0,drop=FALSE])
cosine <- function(a,b) sum(a*b)/sqrt(sum(a^2)*sum(b^2))
cv     <- function(b,l){ t<-table(b,l); sqrt(as.numeric(suppressWarnings(chisq.test(t)$statistic))/(sum(t)*(min(dim(t))-1))) }
two_step <- function(dat,bat,lab,tst,ltst,sup1,sup2){
  dat_c<-if(length(unique(bat))>=2) suppressMessages(ComBat(dat,batch=bat,mod=if(sup1)model.matrix(~lab)else NULL)) else dat
  comb<-cbind(dat_c,tst); cb<-c(rep(1L,ncol(dat_c)),rep(2L,ncol(tst)))
  out<-suppressMessages(ComBat(comb,batch=cb,mod=if(sup2)model.matrix(~c(lab,ltst))else NULL,ref.batch=1L))
  list(ref=dat_c, tgt=out[,(ncol(dat_c)+1):ncol(out),drop=FALSE])
}
cen_bacc <- function(ref,tgt,rl,tl){
  ca<-rowMeans(ref[,rl==1,drop=FALSE]); ch<-rowMeans(ref[,rl==0,drop=FALSE])
  da<-sqrt(colSums((tgt-ca)^2)); dh<-sqrt(colSums((tgt-ch)^2)); pred<-as.integer(da<dh)
  tp<-sum(pred==1&tl==1);tn<-sum(pred==0&tl==0);fp<-sum(pred==1&tl==0);fn<-sum(pred==0&tl==1)
  mean(c(if((tp+fn)>0)tp/(tp+fn)else NA, if((tn+fp)>0)tn/(tn+fp)else NA),na.rm=TRUE)
}
whiten_test <- function(ref,tgt,rl,tl,K=50,ridge=1e-2){
  K<-min(K,ncol(ref)-2); pca<-prcomp(t(ref),rank.=K); Rp<-pca$x; Tp<-predict(pca,t(tgt))
  m1<-colMeans(Rp[rl==1,,drop=FALSE]); m0<-colMeans(Rp[rl==0,,drop=FALSE]); dmu<-m1-m0
  c1<-scale(Rp[rl==1,,drop=FALSE],m1,FALSE); c0<-scale(Rp[rl==0,,drop=FALSE],m0,FALSE)
  Sw<-(crossprod(c1)+crossprod(c0))/(nrow(Rp)-2); Sw<-Sw+diag(ridge*mean(diag(Sw)),K)
  dp<-function(w){ s<-as.numeric(Tp%*%w)
    sp<-sqrt(((sum(tl==1)-1)*var(s[tl==1])+(sum(tl==0)-1)*var(s[tl==0]))/(length(s)-2))
    (mean(s[tl==1])-mean(s[tl==0]))/max(sp,1e-8) }
  c(raw=dp(dmu/sqrt(sum(dmu^2))), whitened=dp({w<-solve(Sw,dmu); w/sqrt(sum(w^2))}))
}
build <- function(n,test){
  refs<-train_studies(n,test); cg<-Reduce(intersect,lapply(c(refs,test),function(s)rownames(dat_lst[[s]])))
  list(dat=log_safe(do.call(cbind,lapply(refs,function(s)dat_lst[[s]][cg,]))),
       bat=do.call(c,lapply(refs,function(s)rep(s,ncol(dat_lst[[s]])))),
       lab=do.call(c,lapply(refs,function(s)bin(label_lst[[s]]))),
       tst=log_safe(dat_lst[[test]][cg,]), ltst=bin(label_lst[[test]]), refs=refs)
}

for (sc in list(c(5,"GSE37250_SA"), c(3,"USA"))) {
  n<-as.integer(sc[1]); test<-sc[2]; B<-build(n,test)
  cat(sprintf("\n========== %s scenario: n=%d test=%s  (Cramer's V=%.3f) ==========\n",
              ifelse(test=="USA","WORKING","FAILING"), n, test, cv(B$bat,B$lab)))
  marg <- wvec(B$dat,B$lab)
  des  <- cbind(model.matrix(~ -1 + factor(B$bat)), class=B$lab)
  Bh   <- solve(crossprod(des), crossprod(des, t(B$dat))); partial <- Bh["class",]
  fs   <- two_step(B$dat,B$bat,B$lab,B$tst,B$ltst,TRUE,FALSE); w_sup<-wvec(fs$ref,B$lab)
  w_test <- wvec(B$tst,B$ltst)
  cat(sprintf("  #12 cos(marginal,partial)=%+.3f  cos(partial,w_sup)=%+.3f  cos(marginal,w_sup)=%+.3f\n",
              cosine(marg,partial), cosine(partial,w_sup), cosine(marg,w_sup)))
  cat(sprintf("      cos(w_sup,TEST)=%+.3f  cos(marginal,TEST)=%+.3f  | sup centroid bacc=%.3f\n",
              cosine(w_sup,w_test), cosine(marg,w_test), cen_bacc(fs$ref,fs$tgt,B$lab,B$ltst)))
  ws <- whiten_test(fs$ref, fs$tgt, B$lab, B$ltst)
  cat(sprintf("  #15 supervised test d':  raw(KNN)=%+.3f   whitened(LDA)=%+.3f\n", ws["raw"], ws["whitened"]))
}

# ── #14 causal de-confounding on the FAILING scenario ─────────────────────────
cat("\n========== #14 de-confounding (n=5 test=GSE37250_SA) ==========\n")
B<-build(5,"GSE37250_SA"); w_test<-wvec(B$tst,B$ltst); marg<-wvec(B$dat,B$lab)
set.seed(1)
keep<-unlist(lapply(B$refs,function(s){idx<-which(B$bat==s);l<-B$lab[idx];k<-min(sum(l==1),sum(l==0))
  c(sample(idx[l==1],k),sample(idx[l==0],k))}))
for(tag in c("orig","balanced")){
  if(tag=="orig"){d<-B$dat;b<-B$bat;l<-B$lab}else{d<-B$dat[,keep];b<-B$bat[keep];l<-B$lab[keep]}
  fs<-two_step(d,b,l,B$tst,B$ltst,TRUE,FALSE); ws<-wvec(fs$ref,l)
  cat(sprintf("  sup-%-8s V=%.3f n=%d | cos(w_sup,TEST)=%+.3f  centroid bacc=%.3f\n",
              tag, cv(b,l), length(l), cosine(ws,w_test), cen_bacc(fs$ref,fs$tgt,l,B$ltst)))
}
fu<-two_step(B$dat,B$bat,B$lab,B$tst,B$ltst,FALSE,FALSE)
cat(sprintf("  unsup      V=%.3f       | cos(w_unsup,TEST)=%+.3f  centroid bacc=%.3f\n",
            cv(B$bat,B$lab), cosine(wvec(fu$ref,B$lab),w_test), cen_bacc(fu$ref,fu$tgt,B$lab,B$ltst)))
cat("\n-> done\n")
