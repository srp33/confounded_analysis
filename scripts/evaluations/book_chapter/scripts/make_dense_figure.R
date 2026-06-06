#!/usr/bin/env Rscript
# make_dense_figure.R — single dense 6-panel walkthrough (3 rows x 2 cols).
# Minimal on-panel text; the narrative lives in the caption (printed + saved .txt).
#
#  ROW 1  SETUP & PROBLEM
#   1  Data & task: 6 TB cohorts, leave-one-study-out, ComBat on the rest, classify held-out.
#   2  Problem: only KNN collapses under supervised ComBat (MCC by classifier).
#  ROW 2  ALGORITHM & MECHANISM
#   3  Algorithm: what ComBat does to ONE real gene — shift studies together, scale spreads,
#      re-inject a class effect (supervised). Per-gene, per-batch, positive-affine.
#   4  Mechanism: per-study leftovers cancel with equal weight, but dividing by sqrt(batch
#      scale) un-balances them -> they pile up into a residual that opposes the class axis.
#  ROW 3  CONSEQUENCE & FIX
#   5  Consequence: held-out samples now sit next to the WRONG class -> KNN wrong neighbours.
#   6  Fix: mean.only (drop the per-gene scale) -> right axis, right neighbours.

suppressMessages({ library(sva); library(matrixStats); library(ggplot2); library(gridExtra); library(grid) })
load("data/TB_real_data.RData"); OUT_DIR<-"outputs/diagnostics/hypothesis_tests"
ALL_STUDIES<-c("GSE37250_SA","USA","India","GSE37250_M","Africa","GSE39941_M")
bin<-function(l) as.integer(ifelse(l%in%c("1",1,"Active"),1L,0L))
train_studies<-function(n,test) ALL_STUDIES[ALL_STUDIES!=test][seq_len(n)]
log_safe<-function(m){if(max(m,na.rm=TRUE)>100){mn<-min(m,na.rm=TRUE);if(mn<0)m<-m-mn;m<-log2(m+1)};m}
unitv<-function(v) v/sqrt(sum(v^2))
clab<-function(l) factor(ifelse(l==1,"Active TB","Control"),levels=c("Control","Active TB"))
ccol<-c("Control"="#4e79a7","Active TB"="#e15759")
adjcol<-c("unsupervised ComBat"="#4e79a7","supervised ComBat"="#e15759")

n<-5; test<-"GSE37250_SA"; refs<-train_studies(n,test)
cg<-Reduce(intersect,lapply(c(refs,test),function(s)rownames(dat_lst[[s]])))
dat<-log_safe(do.call(cbind,lapply(refs,function(s)dat_lst[[s]][cg,])))
bat<-as.factor(do.call(c,lapply(refs,function(s)rep(s,ncol(dat_lst[[s]])))))
lab<-do.call(c,lapply(refs,function(s)bin(label_lst[[s]])))
tst<-log_safe(dat_lst[[test]][cg,]); ltst<-bin(label_lst[[test]]); ubat<-levels(bat); N<-ncol(dat); nb<-length(ubat)

# ── shared formula ingredients ───────────────────────────────────────────────
design<-cbind(model.matrix(~ -1+bat), class=lab); B<-solve(crossprod(design),crossprod(design,t(dat)))
class_effect<-B["class",]; nbi<-as.numeric(table(bat)); grand<-as.numeric(crossprod(nbi/N,B[1:nb,]))
var_pooled<-rowMeans((dat-t(design%*%B))^2)
s_data<-(dat-(outer(grand,rep(1,N))+outer(class_effect,lab)))/sqrt(var_pooled)
delta<-sapply(ubat,function(b) rowVars(s_data[,bat==b,drop=FALSE]))
gap_i<-sapply(ubat,function(b){i<-bat==b; rowMeans(dat[,i&lab==1,drop=FALSE])-rowMeans(dat[,i&lab==0,drop=FALSE])})
n1<-sapply(ubat,function(b)sum(bat==b&lab==1)); n0<-sapply(ubat,function(b)sum(bat==b&lab==0))
p<-n1/(n1+n0); w_i<-(n1/sum(lab==1))*(1-p)+(n0/sum(lab==0))*p
interaction<-sweep(gap_i,1,class_effect,"-")
ref_sup<-suppressMessages(ComBat(dat,batch=bat,mod=model.matrix(~lab)))
ref_mo <-suppressMessages(ComBat(dat,batch=bat,mod=model.matrix(~lab),mean.only=TRUE))
w_sup<-rowMeans(ref_sup[,lab==1])-rowMeans(ref_sup[,lab==0]); residual<-w_sup-class_effect
e1<-unitv(class_effect); e2<-unitv(residual-sum(residual*e1)*e1)

# ════ PANEL 1 — data & task schematic ════════════════════════════════════════
cnt<-do.call(rbind,lapply(ALL_STUDIES,function(s){l<-bin(label_lst[[s]])
  rbind(data.frame(study=s,cls="Control",  nv=sum(l==0)),
        data.frame(study=s,cls="Active TB",nv=sum(l==1)))}))
cnt$study<-factor(cnt$study,levels=rev(ALL_STUDIES)); cnt$cls<-factor(cnt$cls,levels=c("Control","Active TB"))
tot<-tapply(cnt$nv,cnt$study,sum); heldidx<-which(rev(ALL_STUDIES)==test)
ytxt_col<-ifelse(rev(ALL_STUDIES)==test,"#b3261e","grey30")
ytxt_face<-ifelse(rev(ALL_STUDIES)==test,"bold","plain")
p1<-ggplot(cnt,aes(study,nv,fill=cls))+
  geom_col(width=0.72)+
  annotate("segment",x=heldidx,xend=heldidx,y=-Inf,yend=Inf,colour="#b3261e",linewidth=8,alpha=0.10)+
  coord_flip()+
  geom_text(data=data.frame(study=names(tot),nv=tot),aes(study,nv,label=paste0(" n=",nv)),
            inherit.aes=FALSE,hjust=0,size=2.9,colour="grey30")+
  scale_fill_manual(values=ccol,name=NULL)+
  scale_y_continuous(expand=expansion(mult=c(0.02,0.12)))+
  labs(title="1  Six TB blood cohorts; leave one study out",
       subtitle="Batch-correct the other studies together (ComBat), then classify the held-out study (red).",
       x=NULL,y="samples")+
  theme_minimal(base_size=10)+theme(legend.position="top",legend.key.size=unit(0.35,"cm"),
       panel.grid.major.y=element_blank(),
       axis.text.y=element_text(colour=ytxt_col,face=ytxt_face))

# ════ PANEL 2 — problem: MCC by classifier (precomputed) ══════════════════════
bc<-read.csv(file.path(OUT_DIR,"knn_vs_rda_by_classifier.csv"))
m<-bc[bc$metric=="mcc",]
nice<-c(knn="KNN",shrinkageLDA="LDA",svm="SVM",elasticnet="elastic-net",rf="random forest",
        xgboost="XGBoost",nnet="neural net")
m$cl<-factor(nice[m$classifier],levels=nice[order(m$delta)])   # most-affected (KNN) first
ml<-rbind(data.frame(cl=m$cl,adj="unsupervised ComBat",mcc=m$combat),
          data.frame(cl=m$cl,adj="supervised ComBat",  mcc=m$combat_sup))
ml$adj<-factor(ml$adj,levels=names(adjcol))
p2<-ggplot(ml,aes(cl,mcc,fill=adj))+
  geom_hline(yintercept=0,colour="grey55",linewidth=0.4)+
  geom_col(position=position_dodge(width=0.72),width=0.66)+
  scale_fill_manual(values=adjcol,name=NULL)+
  labs(title="2  Only KNN degrades under supervised ComBat",
       subtitle="Held-out classification quality (MCC). Whitening, margin and tree models are unaffected;\nKNN alone falls from 0.63 to below chance (-0.19).",
       x=NULL,y="MCC (held-out)")+
  theme_minimal(base_size=10)+theme(legend.position="top",legend.key.size=unit(0.35,"cm"),
       axis.text.x=element_text(angle=25,hjust=1))

# ════ PANEL 3 — algorithm on ONE real gene ═══════════════════════════════════
flip<-sign(w_sup)!=sign(class_effect) & abs(class_effect)>quantile(abs(class_effect),0.6)
gi<-which(flip)[which.max(abs(class_effect[flip]))]               # a representative class gene that flips
mk<-function(mat,stage) do.call(rbind,lapply(seq_len(nb),function(k){i<-bat==ubat[k]
  data.frame(study=sub("_$","",ubat[k]),expr=mat[gi,i],cls=clab(lab[i]),stage=stage)}))
g3<-rbind(mk(dat,"raw"),mk(ref_sup,"supervised\nComBat"))
g3$stage<-factor(g3$stage,levels=c("raw","supervised\nComBat"))
g3$study<-factor(g3$study,levels=sub("_$","",ubat))
# mark the lowest-scale study (its spread is divided by the smallest sqrt-delta -> amplified)
sd_raw<-sapply(seq_len(nb),function(k) sd(dat[gi,bat==ubat[k]])); lo<-sub("_$","",ubat)[which.min(sd_raw)]
p3plot<-ggplot(g3,aes(study,expr,fill=cls))+
  geom_boxplot(outlier.size=0.4,linewidth=0.3,position=position_dodge(width=0.78),width=0.7)+
  facet_wrap(~stage,nrow=1)+scale_fill_manual(values=ccol,name=NULL)+
  labs(title="3  What ComBat does to one gene",
       subtitle=sprintf("Within each study it shifts means together and rescales spread, then (supervised) re-injects\na class effect. The lowest-variance study (%s) has its rescaled signal amplified and the gap flips.",lo),
       x=NULL,y="expression (log2)")+
  theme_minimal(base_size=10)+theme(legend.position="top",legend.key.size=unit(0.35,"cm"),
       axis.text.x=element_text(angle=25,hjust=1),strip.text=element_text(face="bold"))
eq3<-arrangeGrob(
  textGrob("per-gene correction (gene g, study i):",gp=gpar(cex=0.72,col="grey35",fontface="italic"),x=unit(0.03,"npc"),hjust=0),
  textGrob(expression(gap[g]^{corr}==beta[g]+sum(w[i]*(gap[g*i]-beta[g])/sqrt(delta[g*i]),i)),
           gp=gpar(cex=0.9,col="grey15")),
  textGrob(expression(beta[g]~"= one shared pooled effect      "*delta[g*i]~"= within-study scale (small "*delta~"-> amplified)"),
           gp=gpar(cex=0.7,col="grey40")),
  heights=c(0.3,0.42,0.28))
p3<-arrangeGrob(p3plot,eq3,heights=c(0.78,0.22))

# ════ PANEL 4 — mechanism: per-study breakdown of the broken cancellation ════
# Work along the class direction (1D; orthogonal part omitted). Pooled effect is ONE shared value.
pooled<-sum(class_effect*e1)
gap_proj<-sapply(seq_len(nb),function(k) sum(gap_i[,k]*e1))          # each study's class gap, projected
leftover<-gap_proj-pooled                                           # = (gap_i - pooled), per study
eqw<-w_i*leftover                                                   # equal-weighted contribution  (sum ~ 0)
scl<-sapply(seq_len(nb),function(k) w_i[k]*sum(interaction[,k]/sqrt(pmax(delta[,k],1e-8))*e1)) # /sqrt(scale)
rep_sd<-sapply(seq_len(nb),function(k) sqrt(median(delta[,k])))     # representative within-study scale
sdsp<-sub("_$","",ubat); ordp<-order(gap_proj); levp<-sdsp[ordp]

d4a<-data.frame(study=factor(sdsp,levels=levp),gap=gap_proj,sc=rep_sd)
p4a<-ggplot(d4a,aes(gap,study))+
  geom_vline(xintercept=pooled,linetype="dashed",colour="grey45")+
  annotate("text",x=pooled,y=nb+0.6,label="one shared\npooled effect",hjust=0.5,vjust=0.5,size=2.4,colour="grey45")+
  geom_segment(aes(x=pooled,xend=gap,y=study,yend=study),colour="grey75",linewidth=0.6)+
  geom_point(aes(colour=sc),size=2.8)+
  scale_colour_viridis_c(name="within-study\nscale",option="plasma",end=0.9)+
  scale_y_discrete(expand=expansion(add=c(0.6,1.2)))+
  labs(title="4  Each study's class gap vs the shared pooled effect",
       subtitle="Stem length = leftover (gap - pooled). Colour = within-study scale; small-scale studies\n(dark) get divided by a small sqrt-scale below, amplifying their leftover.",
       x="class gap, projected on the class direction",y=NULL)+
  theme_minimal(base_size=10)+theme(legend.position="right",legend.key.size=unit(0.3,"cm"))

TOT<-"net (sum)"
d4b<-rbind(
  data.frame(study=sdsp,type="equal weighting",val=eqw),
  data.frame(study=sdsp,type="divide by sqrt(scale)",val=scl),
  data.frame(study=TOT, type="equal weighting",val=sum(eqw)),
  data.frame(study=TOT, type="divide by sqrt(scale)",val=sum(scl)))
d4b$study<-factor(d4b$study,levels=c(TOT,levp))
d4b$type<-factor(d4b$type,levels=c("equal weighting","divide by sqrt(scale)"))
b4col<-c("equal weighting"="#9bb8d3","divide by sqrt(scale)"="#e15759")
p4b<-ggplot(d4b,aes(val,study,fill=type))+
  geom_hline(yintercept=1.5,linetype="dotted",colour="grey60")+
  geom_vline(xintercept=0,colour="grey55")+
  geom_col(position=position_dodge(width=0.7),width=0.62)+
  annotate("text",x=sum(eqw),y=1.20,label=sprintf("equal: %.1f (cancels)",sum(eqw)),
           hjust=-0.06,size=2.4,colour="#3b5a78")+
  annotate("text",x=sum(scl),y=0.80,label=sprintf("scaled: %.0f (against class)",sum(scl)),
           hjust=-0.03,size=2.4,colour="#b3261e",fontface="bold")+
  scale_fill_manual(values=b4col,name=NULL)+
  scale_x_continuous(expand=expansion(mult=c(0.06,0.10)))+
  labs(title=NULL,
       subtitle="Equal weighting: contributions cancel (net ~ 0). Dividing by sqrt(scale) up-weights the\nlow-scale studies so they stop cancelling -> a net residual opposing the class.",
       x="contribution to the corrected class gap",y=NULL)+
  theme_minimal(base_size=10)+theme(legend.position="top",legend.key.size=unit(0.35,"cm"))
p4<-arrangeGrob(p4a,p4b,heights=c(0.5,0.5))

# ════ PANELS 5 & 6 — KNN's view: consequence & fix ═══════════════════════════
align<-function(ref){comb<-cbind(ref,tst);cb<-c(rep(1L,ncol(ref)),rep(2L,ncol(tst)))
  o<-suppressMessages(ComBat(comb,batch=cb,mod=NULL,ref.batch=1L)); o[,(ncol(ref)+1):ncol(o),drop=FALSE]}
hull<-function(d) d[chull(d$x,d$y),]
panelCD<-function(ref,ttl,sub){
  tgt<-align(ref); hvg<-order(rowVars(ref),decreasing=TRUE)[1:1000]
  Rh<-ref[hvg,,drop=FALSE]; Th<-tgt[hvg,,drop=FALSE]
  pca<-prcomp(t(Rh),rank.=2); R<-pca$x; Tt<-predict(pca,t(Th))
  dfR<-data.frame(x=R[,1],y=R[,2],cls=clab(lab)); dfT<-data.frame(x=Tt[,1],y=Tt[,2],cls=clab(ltst))
  hl<-rbind(cbind(hull(dfR[dfR$cls=="Active TB",]),g="a"),cbind(hull(dfR[dfR$cls=="Control",]),g="c"))
  cen<-rbind(data.frame(x=mean(dfR$x[lab==1]),y=mean(dfR$y[lab==1]),cls="Active TB"),
             data.frame(x=mean(dfR$x[lab==0]),y=mean(dfR$y[lab==0]),cls="Control"))
  set.seed(3); pick<-sample(ncol(Th),16)
  nn<-sapply(pick,function(j) which.min(colSums((Rh-Th[,j])^2)))
  links<-data.frame(x=Tt[pick,1],y=Tt[pick,2],xe=R[nn,1],ye=R[nn,2],
                    correct=ifelse(lab[nn]==ltst[pick],"neighbour correct","neighbour WRONG"))
  pctw<-round(100*mean(links$correct=="neighbour WRONG"))
  ggplot()+
    geom_polygon(data=hl,aes(x,y,group=g,fill=ifelse(g=="a","Active TB","Control")),alpha=0.12)+
    geom_point(data=dfR,aes(x,y,colour=cls),size=0.8,alpha=0.3)+
    geom_segment(data=links,aes(x,y,xend=xe,yend=ye,linetype=correct),colour="grey20",linewidth=0.45)+
    geom_point(data=dfT,aes(x,y,colour=cls),shape=17,size=2.2,alpha=0.95)+
    geom_point(data=cen,aes(x,y,fill=cls),shape=23,size=4.2,stroke=1,colour="black")+
    scale_colour_manual(values=ccol,name="sample",guide=guide_legend(order=1))+
    scale_fill_manual(values=ccol,guide="none")+
    scale_linetype_manual(values=c("neighbour correct"="solid","neighbour WRONG"="22"),name=NULL)+
    annotate("text",x=Inf,y=Inf,hjust=1.03,vjust=1.4,size=3,fontface="bold",
             colour=ifelse(pctw>50,"#b3261e","#1b7837"),
             label=sprintf("%d%% neighbours WRONG class",pctw))+
    labs(title=ttl,subtitle=sub,x="PC1 of corrected data (what KNN weights most)",y="PC2")+
    theme_minimal(base_size=10)+theme(legend.position="top",legend.key.size=unit(0.35,"cm"))
}
p5<-panelCD(ref_sup,"5  Supervised: held-out samples land by the WRONG class",
            "KNN's view. Triangles = held-out study; diamonds = training class centres.\nActive triangles fall in the training Control cloud (dashed links = wrong neighbour).")
p6<-panelCD(ref_mo,"6  Fix (mean.only): right class, right neighbours",
            "Same study after dropping ONLY the per-gene scale step. Active triangles now\nneighbour training Active (solid links).")

g<-arrangeGrob(p1,p2,p3,p4,p5,p6,ncol=2,
   layout_matrix=matrix(1:6,ncol=2,byrow=TRUE))
ggsave(file.path(OUT_DIR,"dense_figure.png"),g,width=15,height=16.5,dpi=200)
ggsave(file.path(OUT_DIR,"dense_figure.pdf"),g,width=15,height=16.5)
cat("-> dense_figure.png/.pdf\n")
