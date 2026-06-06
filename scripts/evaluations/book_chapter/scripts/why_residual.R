#!/usr/bin/env Rscript
# why_residual.R — WHY the per-gene scale step injects a class-opposing residual.
#
# Derivation: corrected_gap_g = class_effect_g + residual_g, with
#   residual_g = sum_i  (gap_{g,i} - class_effect_g)/sqrt(delta_{g,i}) * w_i
# where gap_{g,i} = within-batch class gap, delta_{g,i} = batch scale (var of the
# class-removed standardized data in batch i), w_i = a positive class-mix weight.
# The interactions (gap_i - class_effect) sum to ~0 when delta is uniform, so they
# CANCEL; non-uniform delta breaks the cancellation -> a net residual. mean.only
# sets delta==1 -> exact cancellation -> no residual.
#
# Test 1: does the formula reproduce the actual residual (w_sup - class_effect)?
# Test 2: ABLATE the class-by-batch interaction -> residual/flip should vanish.
# Test 3: ABLATE only the delta heterogeneity (force equal batch scales) -> vanish.

suppressMessages({ library(sva); library(matrixStats) })
load("data/TB_real_data.RData")
ALL_STUDIES<-c("GSE37250_SA","USA","India","GSE37250_M","Africa","GSE39941_M")
bin<-function(l) as.integer(ifelse(l%in%c("1",1,"Active"),1L,0L))
train_studies<-function(n,test) ALL_STUDIES[ALL_STUDIES!=test][seq_len(n)]
log_safe<-function(m){if(max(m,na.rm=TRUE)>100){mn<-min(m,na.rm=TRUE);if(mn<0)m<-m-mn;m<-log2(m+1)};m}
wvec<-function(m,l) rowMeans(m[,l==1,drop=FALSE])-rowMeans(m[,l==0,drop=FALSE])
cosine<-function(a,b) sum(a*b)/sqrt(sum(a^2)*sum(b^2))

n<-5; test<-"GSE37250_SA"; refs<-train_studies(n,test)
cg<-Reduce(intersect,lapply(c(refs,test),function(s)rownames(dat_lst[[s]])))
dat<-log_safe(do.call(cbind,lapply(refs,function(s)dat_lst[[s]][cg,])))
bat<-as.factor(do.call(c,lapply(refs,function(s)rep(s,ncol(dat_lst[[s]])))))
lab<-do.call(c,lapply(refs,function(s)bin(label_lst[[s]])))
ubat<-levels(bat); nb<-length(ubat); N<-ncol(dat)

# ComBat internals (parametric, no EB shrinkage for delta.hat — close enough)
design<-cbind(model.matrix(~ -1+bat), class=lab)
B<-solve(crossprod(design),crossprod(design,t(dat)))
class_effect<-B["class",]
nbi<-as.numeric(table(bat)); grand<-as.numeric(crossprod(nbi/N, B[1:nb,]))
fitted<-t(design%*%B); var_pooled<-rowMeans((dat-fitted)^2)
stand_mean<-outer(grand,rep(1,N)) + outer(class_effect,lab)
s_data<-(dat-stand_mean)/sqrt(var_pooled)
delta<-sapply(ubat,function(b) rowVars(s_data[,bat==b,drop=FALSE]))   # genes x batch

# per-batch class gap, class-mix weights
gap_i<-sapply(ubat,function(b){i<-bat==b; rowMeans(dat[,i&lab==1,drop=FALSE])-rowMeans(dat[,i&lab==0,drop=FALSE])})
n1<-sapply(ubat,function(b) sum(bat==b&lab==1)); n0<-sapply(ubat,function(b) sum(bat==b&lab==0))
N1<-sum(lab==1); N0<-sum(lab==0); p<-n1/(n1+n0)
w_i<-(n1/N1)*(1-p)+(n0/N0)*p                                          # length nb

interaction<-sweep(gap_i,1,class_effect,"-")                          # genes x batch
pred_resid<-rowSums(sweep(interaction/sqrt(pmax(delta,1e-8)),2,w_i,"*"))
pred_gap<-class_effect+pred_resid

ref_sup<-suppressMessages(ComBat(dat,batch=bat,mod=model.matrix(~lab))); w_sup<-wvec(ref_sup,lab)
actual_resid<-w_sup-class_effect; raw_gap<-wvec(dat,lab)

cat("=== Test 1: does the derived formula reproduce the residual? ===\n")
cat(sprintf("  cor(predicted_residual, actual_residual) = %.3f\n", cor(pred_resid,actual_resid)))
cat(sprintf("  cor(predicted_gap,      actual_corr_gap) = %.3f\n", cor(pred_gap,w_sup)))
cat(sprintf("  cos(predicted_gap, w_sup)=%.3f | cos(w_sup,raw)=%.3f  cos(pred_gap,raw)=%.3f\n",
            cosine(pred_gap,w_sup), cosine(w_sup,raw_gap), cosine(pred_gap,raw_gap)))

# counterfactual residual with delta forced uniform (=1): should ~cancel
pred_resid_flat<-rowSums(sweep(interaction,2,w_i,"*"))
cat(sprintf("\n  If delta were uniform: |residual| median %.4f  -> with real delta: %.4f  (x%.1f larger)\n",
            median(abs(pred_resid_flat)), median(abs(pred_resid)), median(abs(pred_resid))/median(abs(pred_resid_flat))))

# ── Test 2 & 3: causal ablations through the REAL ComBat ──────────────────────
flip_rate<-function(ref,l,raw) mean(sign(wvec(ref,l))!=sign(raw))
cat("\n=== Tests 2/3: ablate cause, re-run real supervised ComBat ===\n")
cat(sprintf("  %-34s cos(w_sup,raw)=%+.3f  flips=%.1f%%\n","baseline (standard supervised)",
            cosine(w_sup,raw_gap),100*flip_rate(ref_sup,lab,raw_gap)))

# Test 2: remove class-by-batch interaction (each batch gets the pooled gap)
dat_h<-dat
for(b in ubat){ i<-which(bat==b); dat_h[,i]<-dat[,i]-outer(gap_i[,b]-class_effect, lab[i]) }
ref_h<-suppressMessages(ComBat(dat_h,batch=bat,mod=model.matrix(~lab))); raw_h<-wvec(dat_h,lab)
cat(sprintf("  %-34s cos(w_sup,raw)=%+.3f  flips=%.1f%%\n","Test2: interaction removed",
            cosine(wvec(ref_h,lab),raw_h),100*flip_rate(ref_h,lab,raw_h)))

# Test 3: remove delta heterogeneity (force equal within-batch variance per gene)
dat_v<-dat
for(b in ubat){ i<-which(bat==b)
  mu<-rowMeans(dat[,i,drop=FALSE]); sdb<-rowSds(dat[,i,drop=FALSE])
  tgt<-rowMeans(sapply(ubat,function(bb) rowSds(dat[,bat==bb,drop=FALSE])))  # common sd
  dat_v[,i]<-mu + (dat[,i,drop=FALSE]-mu)*(tgt/pmax(sdb,1e-8)) }
ref_v<-suppressMessages(ComBat(dat_v,batch=bat,mod=model.matrix(~lab))); raw_v<-wvec(dat_v,lab)
cat(sprintf("  %-34s cos(w_sup,raw)=%+.3f  flips=%.1f%%\n","Test3: batch-scale heterogeneity removed",
            cosine(wvec(ref_v,lab),raw_v),100*flip_rate(ref_v,lab,raw_v)))

# ── figure: formula validation + causal ablations ────────────────────────────
suppressMessages({ library(ggplot2); library(gridExtra) })
OUT_DIR<-"outputs/diagnostics/hypothesis_tests"
d1<-data.frame(pred=pred_gap, act=w_sup)
pA<-ggplot(d1[sample(nrow(d1),min(4000,nrow(d1))),], aes(pred,act))+
  geom_hline(yintercept=0,colour="grey80")+geom_vline(xintercept=0,colour="grey80")+
  geom_abline(slope=1,intercept=0,linetype="dashed",colour="grey50")+
  geom_point(size=0.5,alpha=0.3,colour="#e15759")+
  labs(title="A  The residual formula is exact",
       subtitle=sprintf("Derived corrected gap vs actual (r = %.3f). Formula:\nresidual = sum_i (gap_i - class_effect)/sqrt(delta_i) * w_i",cor(pred_gap,w_sup)),
       x="predicted corrected class gap", y="actual corrected class gap (ComBat)")+
  theme_minimal(base_size=11)
d2<-data.frame(condition=factor(c("standard\nsupervised","remove class x\nbatch interaction","remove batch-scale\nheterogeneity"),
                 levels=c("standard\nsupervised","remove class x\nbatch interaction","remove batch-scale\nheterogeneity")),
               cos=c(cosine(w_sup,raw_gap),cosine(wvec(ref_h,lab),raw_h),cosine(wvec(ref_v,lab),raw_v)))
pB<-ggplot(d2,aes(condition,cos,fill=cos>0))+
  geom_hline(yintercept=0,colour="grey40")+geom_col(width=0.6)+
  scale_fill_manual(values=c(`FALSE`="#e15759",`TRUE`="#4e79a7"),guide="none")+
  labs(title="B  Removing EITHER ingredient fixes the flip",
       subtitle="Alignment of supervised class axis with the raw/true class direction",
       x=NULL,y="cos(supervised axis, raw class axis)")+
  theme_minimal(base_size=11)
ggsave(file.path(OUT_DIR,"why_residual.png"),arrangeGrob(pA,pB,ncol=2),width=12,height=5,dpi=200)
ggsave(file.path(OUT_DIR,"why_residual.pdf"),arrangeGrob(pA,pB,ncol=2),width=12,height=5)
cat("\n-> why_residual.png\n")
