#!/usr/bin/env Rscript
# verify_mg.R
#  PART 1 (focal): what does "flip" decompose into?  corrected_gap = class_effect
#    (re-injected partial coef) + residual_scaling_part. Show which dominates, and
#    whether flipped genes have large class-by-batch interaction.
#  PART 2 (24 scenarios, production KNN): verify combat_sup_mg (sup step1 FULL +
#    mean-only step2 + global scale) vs combat, combat_sup, and our step1-mean.only
#    fix. Prediction: combat_sup_mg only partially helps (step1 still has the scale).

suppressMessages({ library(sva); library(matrixStats); library(class) })
load("data/TB_real_data.RData"); OUT_DIR<-"outputs/diagnostics/hypothesis_tests"
ALL_STUDIES<-c("GSE37250_SA","USA","India","GSE37250_M","Africa","GSE39941_M")
bin<-function(l) as.integer(ifelse(l%in%c("1",1,"Active"),1L,0L))
train_studies<-function(n,test) ALL_STUDIES[ALL_STUDIES!=test][seq_len(n)]
log_safe<-function(m){if(max(m,na.rm=TRUE)>100){mn<-min(m,na.rm=TRUE);if(mn<0)m<-m-mn;m<-log2(m+1)};m}
wvec<-function(m,l) rowMeans(m[,l==1,drop=FALSE])-rowMeans(m[,l==0,drop=FALSE])
cosine<-function(a,b) sum(a*b)/sqrt(sum(a^2)*sum(b^2))

# ── PART 1: decompose the flip (focal failing scenario) ───────────────────────
n<-5; test<-"GSE37250_SA"; refs<-train_studies(n,test)
cg<-Reduce(intersect,lapply(c(refs,test),function(s)rownames(dat_lst[[s]])))
dat<-log_safe(do.call(cbind,lapply(refs,function(s)dat_lst[[s]][cg,])))
bat<-do.call(c,lapply(refs,function(s)rep(s,ncol(dat_lst[[s]])))); lab<-do.call(c,lapply(refs,function(s)bin(label_lst[[s]])))
raw_gap<-wvec(dat,lab)
des<-cbind(model.matrix(~ -1+factor(bat)), class=lab); partial<-solve(crossprod(des),crossprod(des,t(dat)))["class",]
ref_sup<-suppressMessages(ComBat(dat,batch=bat,mod=model.matrix(~lab))); w_sup<-wvec(ref_sup,lab)
residual<-w_sup-partial                            # the non-class-effect part of the corrected gap
# per-gene class-by-batch interaction: SD across batches of the within-batch class gap
perbatch_gap<-sapply(unique(bat),function(b){i<-bat==b; rowMeans(dat[,i&lab==1,drop=FALSE])-rowMeans(dat[,i&lab==0,drop=FALSE])})
interaction_sd<-rowSds(perbatch_gap); maineff<-abs(rowMeans(perbatch_gap))
flip<-sign(w_sup)!=sign(raw_gap)
cat("PART 1 — decomposing the flip (n=5 test=GSE37250_SA):\n")
cat(sprintf("  corrected_gap = class_effect(partial) + residual.  median|partial|=%.3f  median|residual|=%.3f\n",
            median(abs(partial)), median(abs(residual))))
cat(sprintf("  cos(residual, partial)=%+.3f  cos(residual, w_sup)=%+.3f   (residual dominates & opposes class effect)\n",
            cosine(residual,partial), cosine(residual,w_sup)))
cat(sprintf("  class-by-batch interaction (SD of per-batch gap) / |main effect|:  flipped %.2f   non-flipped %.2f\n",
            median(interaction_sd[flip]/pmax(maineff[flip],1e-6)), median(interaction_sd[!flip]/pmax(maineff[!flip],1e-6))))

# ── PART 2: production KNN for combat / combat_sup / combat_sup_mg / step1-fix ─
knn_prod<-function(ref,tgt,rl,tl,n_genes=1000){
  cgi<-intersect(rownames(ref),rownames(tgt)); ref<-ref[cgi,];tgt<-tgt[cgi,]
  hvg<-order(rowVars(ref),decreasing=TRUE)[seq_len(min(n_genes,nrow(ref)))]
  X<-t(ref[hvg,,drop=FALSE]);Y<-t(tgt[hvg,,drop=FALSE]);ns<-nrow(X)
  ks<-c(3,5,7,9,11);ks<-ks[ks<ns];if(!length(ks))ks<-min(3,ns-1)
  acc<-sapply(ks,function(k){f<-cut(seq_len(ns),5,labels=FALSE);ok<-0;tot<-0
    for(fo in 1:5){te<-which(f==fo);tr<-which(f!=fo);if(length(te)&&length(tr)){
      p<-class::knn(X[tr,,drop=FALSE],X[te,,drop=FALSE],rl[tr],k=k);ok<-ok+sum(p==rl[te]);tot<-tot+length(te)}}
    if(tot)ok/tot else 0}); k<-ks[which.max(acc)]
  pred<-as.integer(as.character(class::knn(X,Y,rl,k=k)))
  tp<-sum(pred==1&tl==1);tn<-sum(pred==0&tl==0);fp<-sum(pred==1&tl==0);fn<-sum(pred==0&tl==1)
  den<-sqrt(as.double(tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
  list(bacc=mean(c(if((tp+fn)>0)tp/(tp+fn)else NA,if((tn+fp)>0)tn/(tn+fp)else NA),na.rm=TRUE),mcc=if(den>0)(tp*tn-fp*fn)/den else 0)
}
methods<-c("combat","combat_sup","combat_sup_mg","fix_step1_meanonly")
build_method<-function(m,dat,bat,lab,tst,ltst){
  if(m=="combat"){ rc<-suppressMessages(ComBat(dat,batch=bat,mod=NULL))
    o<-suppressMessages(ComBat(cbind(rc,tst),batch=c(rep(1L,ncol(rc)),rep(2L,ncol(tst))),mod=NULL,ref.batch=1L))
    return(list(ref=rc,tgt=o[,(ncol(rc)+1):ncol(o),drop=FALSE])) }
  if(m=="combat_sup"){ rc<-suppressMessages(ComBat(dat,batch=bat,mod=model.matrix(~lab)))
    o<-suppressMessages(ComBat(cbind(rc,tst),batch=c(rep(1L,ncol(rc)),rep(2L,ncol(tst))),mod=model.matrix(~c(lab,ltst)),ref.batch=1L))
    return(list(ref=rc,tgt=o[,(ncol(rc)+1):ncol(o),drop=FALSE])) }
  if(m=="combat_sup_mg"){ rc<-suppressMessages(ComBat(dat,batch=bat,mod=model.matrix(~lab)))  # sup step1 FULL
    s2<-suppressMessages(ComBat(cbind(rc,tst),batch=c(rep(1L,ncol(rc)),rep(2L,ncol(tst))),mod=NULL,ref.batch=1L,mean.only=TRUE))
    t2<-s2[,(ncol(rc)+1):ncol(s2),drop=FALSE]; rsd<-sd(as.vector(rc)); tsd<-sd(as.vector(t2))
    return(list(ref=rc,tgt=if(rsd>1e-10&&tsd>1e-10) t2*(rsd/tsd) else t2)) }
  if(m=="fix_step1_meanonly"){ rc<-suppressMessages(ComBat(dat,batch=bat,mod=model.matrix(~lab),mean.only=TRUE))
    o<-suppressMessages(ComBat(cbind(rc,tst),batch=c(rep(1L,ncol(rc)),rep(2L,ncol(tst))),mod=NULL,ref.batch=1L))
    return(list(ref=rc,tgt=o[,(ncol(rc)+1):ncol(o),drop=FALSE])) }
}
rows<-list()
for(n in 2:5) for(test in ALL_STUDIES){
  refs<-train_studies(n,test); cgi<-Reduce(intersect,lapply(c(refs,test),function(s)rownames(dat_lst[[s]])))
  d<-log_safe(do.call(cbind,lapply(refs,function(s)dat_lst[[s]][cgi,])))
  b<-do.call(c,lapply(refs,function(s)rep(s,ncol(dat_lst[[s]])))); l<-do.call(c,lapply(refs,function(s)bin(label_lst[[s]])))
  tt<-log_safe(dat_lst[[test]][cgi,]); lt<-bin(label_lst[[test]])
  for(m in methods){ fit<-tryCatch(build_method(m,d,b,l,tt,lt),error=function(e)NULL); if(is.null(fit))next
    kp<-knn_prod(fit$ref,fit$tgt,l,lt); rows[[length(rows)+1]]<-data.frame(n=n,test=test,method=m,bacc=kp$bacc,mcc=kp$mcc) }
  cat(sprintf("  done n=%d test=%s\n",n,test))
}
R<-do.call(rbind,rows); write.csv(R,file.path(OUT_DIR,"verify_mg.csv"),row.names=FALSE)
cat("\nPART 2 — production KNN medians (24 scenarios):\n")
agg<-do.call(rbind,lapply(methods,function(m){s<-R[R$method==m,]
  data.frame(method=m,bacc=round(median(s$bacc),3),mcc=round(median(s$mcc),3),n_below_chance=sum(s$mcc<0))}))
print(agg,row.names=FALSE); cat("\n-> verify_mg.csv\n")
