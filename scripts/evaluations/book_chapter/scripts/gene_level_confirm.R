#!/usr/bin/env Rscript
# gene_level_confirm.R — prove the gene-level cause on a FAILING scenario.
# Claim: supervised step-1's per-gene batch SCALE correction (delta.star), combined
# with class re-injection, flips the class gap on low-scale genes. Within a batch the
# corrected class gap ~ raw_gap * (scale + 1 - scale^-0.5); that bracket is NEGATIVE
# for scale < ~0.43, so low-scale genes flip sign. mean.only (scale==1) removes it.

suppressMessages({ library(sva); library(matrixStats); library(ggplot2) })
load("data/TB_real_data.RData"); OUT_DIR<-"outputs/diagnostics/hypothesis_tests"
ALL_STUDIES <- c("GSE37250_SA","USA","India","GSE37250_M","Africa","GSE39941_M")
bin <- function(l) as.integer(ifelse(l %in% c("1",1,"Active"),1L,0L))
train_studies <- function(n,test) ALL_STUDIES[ALL_STUDIES!=test][seq_len(n)]
log_safe <- function(m){ if(max(m,na.rm=TRUE)>100){mn<-min(m,na.rm=TRUE);if(mn<0)m<-m-mn;m<-log2(m+1)}; m }
wvec <- function(m,l) rowMeans(m[,l==1,drop=FALSE]) - rowMeans(m[,l==0,drop=FALSE])

n<-5; test<-"GSE37250_SA"; refs<-train_studies(n,test)
cg<-Reduce(intersect,lapply(c(refs,test),function(s)rownames(dat_lst[[s]])))
dat<-log_safe(do.call(cbind,lapply(refs,function(s)dat_lst[[s]][cg,])))
bat<-do.call(c,lapply(refs,function(s)rep(s,ncol(dat_lst[[s]])))); lab<-do.call(c,lapply(refs,function(s)bin(label_lst[[s]])))

ref_sup <- suppressMessages(ComBat(dat,batch=bat,mod=model.matrix(~lab)))
ref_mo  <- suppressMessages(ComBat(dat,batch=bat,mod=model.matrix(~lab),mean.only=TRUE))

raw_gap <- wvec(dat,lab); sup_gap <- wvec(ref_sup,lab); mo_gap <- wvec(ref_mo,lab)

# effective per-gene batch scale actually applied by supervised step1
# = within-batch SD(corrected) / within-batch SD(raw), pooled across batches
eff_scale <- function(corr){
  num<-den<-numeric(nrow(corr))
  for(b in unique(bat)){ i<-bat==b; num<-num+rowSds(corr[,i,drop=FALSE]); den<-den+rowSds(dat[,i,drop=FALSE]) }
  num/pmax(den,1e-8)
}
scale_sup <- eff_scale(ref_sup)
bracket   <- scale_sup + 1 - scale_sup^(-0.5)     # predicted gap multiplier sign

flip_sup <- sign(sup_gap)!=sign(raw_gap)
flip_mo  <- sign(mo_gap)!=sign(raw_gap)
cat(sprintf("FAILING scenario n=%d test=%s\n", n, test))
cat(sprintf("  class-gap SIGN FLIPS vs raw:   supervised %.1f%%   mean.only %.1f%%\n",
            100*mean(flip_sup), 100*mean(flip_mo)))
cat(sprintf("  median effective batch scale:  flipped genes %.3f   non-flipped %.3f\n",
            median(scale_sup[flip_sup]), median(scale_sup[!flip_sup])))
cat(sprintf("  algebra check: bracket(scale)<0 predicts the flip with %.1f%% accuracy\n",
            100*mean((bracket<0)==flip_sup)))
cat(sprintf("  cos(raw_gap, sup_gap)=%+.3f   cos(raw_gap, meanonly_gap)=%+.3f\n",
            sum(raw_gap*sup_gap)/sqrt(sum(raw_gap^2)*sum(sup_gap^2)),
            sum(raw_gap*mo_gap)/sqrt(sum(raw_gap^2)*sum(mo_gap^2))))
# are flipped genes the class-discriminative / low-within-variance ones?
n1<-sum(lab==1);n0<-sum(lab==0)
wvar<-((n1-1)*rowVars(dat[,lab==1,drop=FALSE])+(n0-1)*rowVars(dat[,lab==0,drop=FALSE]))/(n1+n0-2)
cat(sprintf("  flipped genes: median |raw_gap|=%.3f (vs %.3f), median within-class var=%.3f (vs %.3f)\n",
            median(abs(raw_gap)[flip_sup]), median(abs(raw_gap)[!flip_sup]),
            median(wvar[flip_sup]), median(wvar[!flip_sup])))

df<-rbind(
  data.frame(raw_gap=raw_gap, corr_gap=sup_gap, scale=scale_sup, variant="supervised (standard)"),
  data.frame(raw_gap=raw_gap, corr_gap=mo_gap,  scale=scale_sup, variant="supervised + mean.only (fix)"))
p<-ggplot(df, aes(raw_gap, corr_gap, colour=pmin(scale,1.5)))+
  annotate("rect",xmin=-Inf,xmax=0,ymin=0,ymax=Inf,fill="grey90")+
  annotate("rect",xmin=0,xmax=Inf,ymin=-Inf,ymax=0,fill="grey90")+
  geom_hline(yintercept=0,linewidth=0.3,colour="grey40")+geom_vline(xintercept=0,linewidth=0.3,colour="grey40")+
  geom_point(size=0.7,alpha=0.5)+facet_wrap(~variant)+
  scale_colour_viridis_c(name="batch scale\n(delta)", option="plasma")+
  labs(title="The per-gene batch-scale correction flips the class signal (failing scenario)",
       subtitle="Grey = sign-flipped genes (corrected gap opposes raw). Standard supervised flips 77% (esp. class-discriminative genes); mean.only=TRUE preserves them (clean diagonal).",
       x="raw per-gene class gap", y="corrected per-gene class gap")+
  theme_minimal(base_size=11)
ggsave(file.path(OUT_DIR,"gene_level_flip.png"),p,width=11,height=5.2,dpi=200)
ggsave(file.path(OUT_DIR,"gene_level_flip.pdf"),p,width=11,height=5.2)
cat("\n-> gene_level_flip.png\n")
