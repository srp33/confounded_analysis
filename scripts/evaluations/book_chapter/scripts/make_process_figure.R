#!/usr/bin/env Rscript
# make_process_figure.R — a visual, data-driven walkthrough of the mechanism.
# 4 panels, each a 2D plot where BOTH axes carry information:
#  A  per-study class directions DISAGREE (class x batch interaction)   [raw PCA]
#  B  ComBat's per-gene SCALE breaks the cancellation -> a residual that
#     opposes the class direction                                       [vector sum]
#  C  corrected clouds separate on the WRONG axis; a held-out study's samples
#     sit next to the wrong class -> KNN picks wrong neighbours          [sup PCA]
#  D  mean.only removes the scale step -> right axis, right neighbours   [fix PCA]

suppressMessages({ library(sva); library(matrixStats); library(ggplot2); library(gridExtra) })
load("data/TB_real_data.RData"); OUT_DIR<-"outputs/diagnostics/hypothesis_tests"
ALL_STUDIES<-c("GSE37250_SA","USA","India","GSE37250_M","Africa","GSE39941_M")
bin<-function(l) as.integer(ifelse(l%in%c("1",1,"Active"),1L,0L))
train_studies<-function(n,test) ALL_STUDIES[ALL_STUDIES!=test][seq_len(n)]
log_safe<-function(m){if(max(m,na.rm=TRUE)>100){mn<-min(m,na.rm=TRUE);if(mn<0)m<-m-mn;m<-log2(m+1)};m}
unitv<-function(v) v/sqrt(sum(v^2))

n<-5; test<-"GSE37250_SA"; refs<-train_studies(n,test)
cg<-Reduce(intersect,lapply(c(refs,test),function(s)rownames(dat_lst[[s]])))
dat<-log_safe(do.call(cbind,lapply(refs,function(s)dat_lst[[s]][cg,])))
bat<-as.factor(do.call(c,lapply(refs,function(s)rep(s,ncol(dat_lst[[s]])))))
lab<-do.call(c,lapply(refs,function(s)bin(label_lst[[s]])))
tst<-log_safe(dat_lst[[test]][cg,]); ltst<-bin(label_lst[[test]]); ubat<-levels(bat); N<-ncol(dat); nb<-length(ubat)
clab<-function(l) factor(ifelse(l==1,"Active TB","Control"),levels=c("Control","Active TB"))
ccol<-c("Control"="#4e79a7","Active TB"="#e15759")

# ── shared quantities (the formula ingredients) ──────────────────────────────
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
w_sup<-rowMeans(ref_sup[,lab==1])-rowMeans(ref_sup[,lab==0])
residual<-w_sup-class_effect
e1<-unitv(class_effect); e2<-unitv(residual-sum(residual*e1)*e1)        # plane: class dir x residual dir

# ════ Panels A & B: identical head-to-tail construction; only the weighting differs ════
# leftover_i = gap_i - pooled class effect (= class x batch interaction). Add them
# head-to-tail in the (class direction, orthogonal) plane. A: equal weights -> back
# to ~0 (cancel). B: each divided by sqrt(batch scale) -> pile up into the residual.
proj<-function(v) c(x=sum(v*e1),y=sum(v*e2))
ord<-order(-w_i); studlab<-sub("_$","",ubat)[ord]
chain<-function(scaled){ P<-matrix(0,nb+1,2); for(k in seq_len(nb)){ i<-ord[k]
    v<-interaction[,i]*w_i[i]; if(scaled) v<-interaction[,i]/sqrt(pmax(delta[,i],1e-8))*w_i[i]
    P[k+1,]<-P[k,]+proj(v)}; P}
Pu<-chain(FALSE); Pc<-chain(TRUE)
mkseg<-function(P) data.frame(x=P[-nrow(P),1],y=P[-nrow(P),2],xe=P[-1,1],ye=P[-1,2],
                              mx=(P[-nrow(P),1]+P[-1,1])/2, my=(P[-nrow(P),2]+P[-1,2])/2, study=studlab)
segU<-mkseg(Pu); segC<-mkseg(Pc)

chainPlot<-function(seg,Pend,col,endcol,ttl,sub,endlab,refarrow){
  rx<-range(c(seg$x,seg$xe,0,if(refarrow)max(abs(c(seg$x,seg$xe)))*0.5 else 0))
  g<-ggplot()+geom_hline(yintercept=0,colour="grey88")+geom_vline(xintercept=0,colour="grey88")
  if(refarrow){ L<-max(abs(c(seg$x,seg$xe)))*0.55
    g<-g+annotate("segment",x=0,xend=L,y=0,yend=0,arrow=arrow(length=unit(0.2,"cm"),type="closed"),
                  colour="grey55",linewidth=0.9)+
      annotate("text",x=L,y=0,label="class direction",hjust=1,vjust=-0.8,size=2.9,colour="grey40") }
  g+geom_segment(data=seg,aes(x,y,xend=xe,yend=ye),colour=col,linewidth=1.1,
               arrow=arrow(length=unit(0.14,"cm"),type="closed"))+
    geom_text(data=seg,aes(mx,my,label=study),size=2.5,colour=col,vjust=-0.6,fontface="italic")+
    annotate("point",x=Pend[1],y=Pend[2],size=3.5,colour=endcol)+
    annotate("text",x=Pend[1],y=Pend[2],label=endlab,hjust=0.5,vjust=ifelse(Pend[2]>=0,-0.5,1.3),
             size=3,fontface="bold",colour=endcol)+
    coord_equal()+labs(title=ttl,subtitle=sub,
      x="component along class direction  (-> with class,  <- against)",y="orthogonal component")+
    theme_minimal(base_size=10)
}
maxU<-max(sqrt(rowSums(Pu^2))); distC<-sqrt(sum(Pc[nb+1,]^2)); zoomf<-max(1,round(distC/max(maxU,1e-6)))
pA<-chainPlot(segU,Pu[nb+1,],"#4e79a7","#2f5597",
  "A  Normally the per-study leftovers cancel",
  sprintf("Each study's leftover (its class gap minus the pooled effect) added head-to-tail,\nweighted equally. The walk returns to the origin: the disagreements cancel.\n(axes zoomed ~%dx vs panel B - the walk barely leaves the origin.)",zoomf),
  "Sum ~ 0\n(cancels)",FALSE)
pB<-chainPlot(segC,Pc[nb+1,],"#e15759","#b3261e",
  "B  The per-gene scale breaks the cancellation",
  "Exactly the same leftovers, but each first divided by sqrt(its batch scale).\nNow they pile up - AGAINST the class direction - into a large residual.",
  "RESIDUAL\n(opposes class -> axis flips)",TRUE)

# ════ Panels C/D: corrected samples in the fixed (true class axis, residual axis) frame ════
# x = projection on TRUE class direction (e1); y = projection on residual direction (e2).
align<-function(ref){comb<-cbind(ref,tst);cb<-c(rep(1L,ncol(ref)),rep(2L,ncol(tst)))
  o<-suppressMessages(ComBat(comb,batch=cb,mod=NULL,ref.batch=1L)); o[,(ncol(ref)+1):ncol(o),drop=FALSE]}
hull<-function(d) d[chull(d$x,d$y),]
# KNN's actual view: top-2 PCs of the 1000-HVG corrected reference; neighbours found
# in the real 1000-gene Euclidean distance (what class::knn uses).
panelCD<-function(ref,ttl,sub){
  tgt<-align(ref); hvg<-order(rowVars(ref),decreasing=TRUE)[1:1000]
  Rh<-ref[hvg,,drop=FALSE]; Th<-tgt[hvg,,drop=FALSE]
  pca<-prcomp(t(Rh),rank.=2); R<-pca$x; Tt<-predict(pca,t(Th))
  if(cor(R[lab==1,1],rep(1,sum(lab==1)))||TRUE){}                       # keep PC sign as-is
  dfR<-data.frame(x=R[,1],y=R[,2],cls=clab(lab)); dfT<-data.frame(x=Tt[,1],y=Tt[,2],cls=clab(ltst))
  hl<-rbind(cbind(hull(dfR[dfR$cls=="Active TB",]),g="a"),cbind(hull(dfR[dfR$cls=="Control",]),g="c"))
  cen<-rbind(data.frame(x=mean(dfR$x[lab==1]),y=mean(dfR$y[lab==1]),cls="Active TB"),
             data.frame(x=mean(dfR$x[lab==0]),y=mean(dfR$y[lab==0]),cls="Control"))
  set.seed(3); pick<-sample(ncol(Th),16)
  nn<-sapply(pick,function(j) which.min(colSums((Rh-Th[,j])^2)))        # nearest in 1000-gene space
  links<-data.frame(x=Tt[pick,1],y=Tt[pick,2],xe=R[nn,1],ye=R[nn,2],
                    correct=ifelse(lab[nn]==ltst[pick],"neighbour correct","neighbour WRONG"))
  pctw<-round(100*mean(links$correct=="neighbour WRONG"))
  ggplot()+
    geom_polygon(data=hl,aes(x,y,group=g,fill=ifelse(g=="a","Active TB","Control")),alpha=0.12)+
    geom_point(data=dfR,aes(x,y,colour=cls),size=0.8,alpha=0.3)+
    geom_segment(data=links,aes(x,y,xend=xe,yend=ye,linetype=correct),colour="grey20",linewidth=0.45)+
    geom_point(data=dfT,aes(x,y,colour=cls),shape=17,size=2.2,alpha=0.95)+
    geom_point(data=cen,aes(x,y,fill=cls),shape=23,size=4.5,stroke=1,colour="black")+
    scale_colour_manual(values=ccol,name="sample",guide=guide_legend(order=1))+
    scale_fill_manual(values=ccol,guide="none")+
    scale_linetype_manual(values=c("neighbour correct"="solid","neighbour WRONG"="22"),name=NULL)+
    annotate("text",x=Inf,y=Inf,hjust=1.03,vjust=1.4,size=3.2,fontface="bold",
             colour=ifelse(pctw>50,"#b3261e","#1b7837"),
             label=sprintf("%d%% of shown neighbours: WRONG class",pctw))+
    labs(title=ttl,subtitle=sub,
         x="PC1 of corrected data (high-variance: what KNN weights)",y="PC2 of corrected data")+
    theme_minimal(base_size=10)+theme(legend.position="right",legend.key.size=unit(0.4,"cm"))
}
ref_mo<-suppressMessages(ComBat(dat,batch=bat,mod=model.matrix(~lab),mean.only=TRUE))
pC<-panelCD(ref_sup,"C  Supervised: held-out samples land by the WRONG class",
            "KNN's view (top variance directions). Triangles = held-out study; filled diamonds =\ntraining class centres. Active triangles sit in the training Control cloud (dashed = wrong).")
pD<-panelCD(ref_mo,"D  Fix (mean.only): held-out samples land by the RIGHT class",
            "Same held-out study after removing only the per-gene scale step: Active triangles\nnow neighbour training Active (solid links).")

g<-arrangeGrob(pA,pB,pC,pD,ncol=2)
ggsave(file.path(OUT_DIR,"process_figure.png"),g,width=15,height=11,dpi=200)
ggsave(file.path(OUT_DIR,"process_figure.pdf"),g,width=15,height=11)
cat("-> process_figure.png/.pdf\n")
