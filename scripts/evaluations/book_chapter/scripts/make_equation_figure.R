#!/usr/bin/env Rscript
# make_equation_figure.R — equation-centric mechanism figure.
# The corrected-class-gap formula sits in the center of the page;
# each surrounding panel illustrates exactly one term of the formula.
#
# Layout (4 rows × 6 units):
#  Row 1: [A: study overview]            [B: classifier performance]
#  Row 2: [C: gap_{g,i}]   [D: beta_g]  [E: delta_{g,i}]
#  Row 3: [F: leftover]  [EQUATION]  [G: amplified leftover]
#  Row 4: [H: sum / cancellation]  [I: KNN fails]  [J: KNN fix]
#
# All gene-level panels (C-G) use gene RPS4Y1 (or the selected focal gene).
# Panel labels A-J are printed in panel titles, color-coded by equation term.

suppressMessages({
  library(sva); library(matrixStats); library(ggplot2)
  library(gridExtra); library(grid)
})
load("data/TB_real_data.RData"); OUT_DIR <- "outputs/diagnostics/hypothesis_tests"
ALL_STUDIES <- c("GSE37250_SA","USA","India","GSE37250_M","Africa","GSE39941_M")
bin       <- function(l) as.integer(ifelse(l %in% c("1",1,"Active"), 1L, 0L))
log_safe  <- function(m) {if(max(m,na.rm=TRUE)>100){mn<-min(m,na.rm=TRUE);if(mn<0)m<-m-mn;m<-log2(m+1)};m}
unitv     <- function(v) v/sqrt(sum(v^2))
clab      <- function(l) factor(ifelse(l==1,"Active TB","Control"), levels=c("Control","Active TB"))
ccol      <- c("Control"="#4e79a7","Active TB"="#e15759")

# ── term colors (used consistently across panels and the equation card) ──────
tc <- c(gap="#2166ac", beta="#b2182b", delta="#1a9641", resid="#762a83")

# ── data setup (n=5, test=GSE37250_SA, the canonical failing scenario) ───────
n <- 5; test <- "GSE37250_SA"
refs <- ALL_STUDIES[ALL_STUDIES != test][seq_len(n)]
cg   <- Reduce(intersect, lapply(c(refs,test), function(s) rownames(dat_lst[[s]])))
dat  <- log_safe(do.call(cbind, lapply(refs, function(s) dat_lst[[s]][cg,])))
bat  <- as.factor(do.call(c, lapply(refs, function(s) rep(s, ncol(dat_lst[[s]])))))
lab  <- do.call(c, lapply(refs, function(s) bin(label_lst[[s]])))
tst  <- log_safe(dat_lst[[test]][cg,]); ltst <- bin(label_lst[[test]])
ubat <- levels(bat); nb <- length(ubat); N <- ncol(dat)

# ── shared ComBat quantities ─────────────────────────────────────────────────
design       <- cbind(model.matrix(~ -1+bat), class=lab)
B            <- solve(crossprod(design), crossprod(design, t(dat)))
class_effect <- B["class",]
nbi          <- as.numeric(table(bat)); grand <- as.numeric(crossprod(nbi/N, B[1:nb,]))
var_pooled   <- rowMeans((dat - t(design %*% B))^2)
stand_mean   <- outer(grand, rep(1,N)) + outer(class_effect, lab)
s_data       <- (dat - stand_mean) / sqrt(var_pooled)
delta        <- sapply(ubat, function(b) rowVars(s_data[, bat==b, drop=FALSE]))  # genes × batches
gap_i        <- sapply(ubat, function(b) {
  i <- bat==b
  rowMeans(dat[,i&lab==1,drop=FALSE]) - rowMeans(dat[,i&lab==0,drop=FALSE])
})
n1_b <- sapply(ubat, function(b) sum(bat==b & lab==1))
n0_b <- sapply(ubat, function(b) sum(bat==b & lab==0))
p_b  <- n1_b/(n1_b+n0_b)
w_i  <- (n1_b/sum(lab==1))*(1-p_b) + (n0_b/sum(lab==0))*p_b
interaction  <- sweep(gap_i, 1, class_effect, "-")

ref_sup <- suppressMessages(ComBat(dat, batch=bat, mod=model.matrix(~lab)))
ref_mo  <- suppressMessages(ComBat(dat, batch=bat, mod=model.matrix(~lab), mean.only=TRUE))
w_sup   <- rowMeans(ref_sup[,lab==1]) - rowMeans(ref_sup[,lab==0])
residual <- w_sup - class_effect
e1 <- unitv(class_effect); e2 <- unitv(residual - sum(residual*e1)*e1)

# ── focal gene selection: large |residual|/|class_effect| AND heterogeneous delta ─
flip      <- sign(w_sup) != sign(class_effect) & abs(class_effect) > quantile(abs(class_effect), 0.6)
delta_cv  <- apply(delta, 1, function(x) sd(x)/pmax(mean(x),1e-8))
score     <- (abs(residual)/pmax(abs(class_effect),1e-6)) * delta_cv
gi        <- which(flip)[which.max(score[flip])]
gene_name <- rownames(dat)[gi]

# per-gene scalar quantities for gene gi
gap_gi      <- gap_i[gi, ]                                  # per-study class gap
beta_gi     <- class_effect[gi]                             # shared pooled estimate
delta_gi    <- delta[gi, ]                                  # per-study scale
leftover_gi <- gap_gi - beta_gi                             # gap_gi - beta_g
amplif_gi   <- leftover_gi / sqrt(pmax(delta_gi, 1e-8))    # ÷ sqrt(delta)

sdsp    <- sub("_$","", ubat)
sg_ord  <- sdsp[order(gap_gi)]                             # study order: ascending gap
sg_fac  <- factor(sdsp, levels=sg_ord)

# ── theme helpers ────────────────────────────────────────────────────────────
th  <- theme_minimal(base_size=9.5) + theme(legend.key.size=unit(0.32,"cm"))
thl <- th + theme(axis.text.y=element_text(size=8))
tag_text <- function(lbl, col) {
  annotate("text", x=-Inf, y=Inf, hjust=-0.25, vjust=1.5,
           label=lbl, size=3.4, fontface="bold", colour=col)
}

# ════════════════════════════════════════════════════════════════════════════
# PANEL A — study overview
# ════════════════════════════════════════════════════════════════════════════
cnt <- do.call(rbind, lapply(ALL_STUDIES, function(s) {
  l <- bin(label_lst[[s]])
  rbind(data.frame(study=s, cls="Control",   nv=sum(l==0)),
        data.frame(study=s, cls="Active TB", nv=sum(l==1)))
}))
cnt$study   <- factor(cnt$study, levels=rev(ALL_STUDIES))
cnt$cls     <- factor(cnt$cls,   levels=c("Control","Active TB"))
tot         <- tapply(cnt$nv, cnt$study, sum)
heldidx     <- which(rev(ALL_STUDIES)==test)
ytxt_col    <- ifelse(rev(ALL_STUDIES)==test, "#b3261e", "grey30")
ytxt_face   <- ifelse(rev(ALL_STUDIES)==test, "bold",    "plain")

pA <- ggplot(cnt, aes(study, nv, fill=cls)) +
  geom_col(width=0.72) +
  annotate("segment", x=heldidx, xend=heldidx, y=-Inf, yend=Inf,
           colour="#b3261e", linewidth=7, alpha=0.09) +
  coord_flip() +
  geom_text(data=data.frame(study=names(tot), nv=tot),
            aes(study, nv, label=paste0(" n=",nv)),
            inherit.aes=FALSE, hjust=0, size=2.7, colour="grey35") +
  scale_fill_manual(values=ccol, name=NULL) +
  scale_y_continuous(expand=expansion(mult=c(0.02,0.12))) +
  labs(title="A  Six TB blood cohorts", x=NULL, y="samples",
       subtitle="Held-out (red): batch-correct the other 5, classify this one.") +
  th + theme(legend.position="top", panel.grid.major.y=element_blank(),
             axis.text.y=element_text(colour=ytxt_col, face=ytxt_face))

# ════════════════════════════════════════════════════════════════════════════
# PANEL B — classifier performance
# ════════════════════════════════════════════════════════════════════════════
bc   <- read.csv(file.path(OUT_DIR,"knn_vs_rda_by_classifier.csv"))
m    <- bc[bc$metric=="mcc",]
nice <- c(knn="KNN",shrinkageLDA="LDA",svm="SVM",elasticnet="elastic-net",
          rf="random forest",xgboost="XGBoost",nnet="neural net")
m$cl <- factor(nice[m$classifier], levels=nice[order(m$delta)])
ml   <- rbind(data.frame(cl=m$cl, adj="unsupervised\nComBat", mcc=m$combat),
              data.frame(cl=m$cl, adj="supervised\nComBat",   mcc=m$combat_sup))
adjcol <- c("unsupervised\nComBat"="#4e79a7","supervised\nComBat"="#e15759")
ml$adj <- factor(ml$adj, levels=names(adjcol))
pB <- ggplot(ml, aes(cl, mcc, fill=adj)) +
  geom_hline(yintercept=0, colour="grey55", linewidth=0.35) +
  geom_col(position=position_dodge(width=0.72), width=0.66) +
  scale_fill_manual(values=adjcol, name=NULL) +
  labs(title="B  Only KNN degrades under supervised ComBat",
       subtitle="MCC averaged across held-out studies. Whitening, margin and tree classifiers\nare unaffected; KNN alone falls from 0.63 to below chance (-0.19).",
       x=NULL, y="MCC") +
  th + theme(legend.position="top", axis.text.x=element_text(angle=25, hjust=1))

# ════════════════════════════════════════════════════════════════════════════
# PANEL C — gap_{g,i}: per-study raw class gap for gene RPS4Y1
# ════════════════════════════════════════════════════════════════════════════
d_c <- data.frame(study=sg_fac, val=gap_gi)
pC <- ggplot(d_c, aes(val, study)) +
  geom_vline(xintercept=0, colour="grey55", linewidth=0.35) +
  geom_segment(aes(x=0, xend=val, y=study, yend=study),
               colour=tc["gap"], linewidth=0.9, alpha=0.7) +
  geom_point(size=2.8, colour=tc["gap"]) +
  geom_text(aes(label=sprintf("%+.2f", val)), hjust=-0.25, size=2.4, colour=tc["gap"]) +
  scale_x_continuous(expand=expansion(mult=c(0.08,0.22))) +
  labs(title=sprintf("C  gap(g,i)  [gene: %s]", gene_name),
       subtitle="Active minus Control mean expression per study.",
       x="class gap (log2)", y=NULL) +
  tag_text("gap(g,i)", tc["gap"]) +
  thl

# ════════════════════════════════════════════════════════════════════════════
# PANEL D — beta_g: the one shared pooled estimate
# ════════════════════════════════════════════════════════════════════════════
d_d <- data.frame(study=sg_fac, gap=gap_gi)
beta_lab <- sprintf("beta = %+.3f", beta_gi)
pD <- ggplot(d_d, aes(gap, study)) +
  geom_vline(xintercept=beta_gi, linetype="dashed",
             colour=tc["beta"], linewidth=0.9) +
  annotate("text", x=beta_gi, y=nb+0.75,
           label=beta_lab, colour=tc["beta"], size=2.8, hjust=0.5, fontface="bold") +
  geom_segment(aes(x=beta_gi, xend=gap, y=study, yend=study),
               colour="grey70", linewidth=0.6) +
  geom_point(size=2.8, colour=tc["gap"]) +
  scale_x_continuous(expand=expansion(mult=c(0.08,0.08))) +
  scale_y_discrete(expand=expansion(add=c(0.5,1.1))) +
  labs(title=sprintf("D  beta_g: one shared pooled estimate  [%s]", gene_name),
       subtitle="beta_g (dashed) = LS coefficient from expr ~ batch + class.\nStem = per-study deviation (the leftover).",
       x="class gap (log2)", y=NULL) +
  tag_text("beta_g", tc["beta"]) +
  thl

# ════════════════════════════════════════════════════════════════════════════
# PANEL E — delta_{g,i}: within-study scale for gene RPS4Y1
# ════════════════════════════════════════════════════════════════════════════
d_e <- data.frame(study=sg_fac, delta=delta_gi, inv_sqrt=1/sqrt(pmax(delta_gi,1e-8)))
pE <- ggplot(d_e, aes(delta, study)) +
  geom_vline(xintercept=1, colour="grey65", linewidth=0.35, linetype="dashed") +
  geom_col(fill=tc["delta"], alpha=0.75, width=0.6) +
  geom_text(aes(label=sprintf("1/sqrt(d)= %+.1f", inv_sqrt)),
            hjust=-0.08, size=2.3, colour=tc["delta"]) +
  scale_x_continuous(expand=expansion(mult=c(0.02,0.30))) +
  labs(title=sprintf("E  delta(g,i): within-study scale  [%s]", gene_name),
       subtitle="Variance of standardized residuals within each study.\nSmall delta -> large 1/sqrt(delta) -> leftover amplified.",
       x="within-study scale (delta)", y=NULL) +
  tag_text("delta(g,i)", tc["delta"]) +
  thl

# ════════════════════════════════════════════════════════════════════════════
# PANEL F — leftover = gap_{g,i} - beta_g (BEFORE dividing by scale)
# ════════════════════════════════════════════════════════════════════════════
d_f <- data.frame(study=sg_fac, val=leftover_gi,
                  wt_contrib=w_i[match(sdsp, sg_ord)] * leftover_gi)
# recompute w_i indexed correctly
wi_named <- setNames(w_i, sdsp)
d_f$wt_contrib <- wi_named[as.character(sg_fac)] * leftover_gi

net_eq <- sum(wi_named * leftover_gi)
pF <- ggplot(d_f, aes(val, study, fill=val>0)) +
  geom_vline(xintercept=0, colour="grey55", linewidth=0.35) +
  geom_col(width=0.65, alpha=0.8) +
  annotate("text", x=0, y=0.4,
           label=sprintf("net sum ~ %.2f (cancels)", net_eq),
           hjust=0.5, vjust=0, size=2.5, colour="grey35", fontface="italic") +
  scale_fill_manual(values=c("TRUE"=tc["gap"], "FALSE"="#d6604d"), guide="none") +
  scale_x_continuous(expand=expansion(mult=c(0.12,0.12))) +
  labs(title=sprintf("F  gap(g,i) - beta_g: leftover  [%s]", gene_name),
       subtitle="Signed deviation of each study's gap from the pooled estimate.\nWith equal weights these roughly cancel (net ~ 0).",
       x="leftover (log2)", y=NULL) +
  tag_text("gap - beta", tc["gap"]) +
  thl

# ════════════════════════════════════════════════════════════════════════════
# PANEL — EQUATION (center)
# ════════════════════════════════════════════════════════════════════════════
mk_card <- function(xmin, xmax, yb, yt, col, sym_lab, desc_lab, panel_lab) {
  list(
    annotate("rect", xmin=xmin, xmax=xmax, ymin=yb, ymax=yt, fill=col, alpha=0.13, colour=col, linewidth=0.4),
    annotate("text", x=(xmin+xmax)/2, y=(yb+yt)*0.5+0.015, label=sym_lab,
             parse=TRUE, size=3.4, colour=col, fontface="bold", hjust=0.5),
    annotate("text", x=(xmin+xmax)/2, y=(yb+yt)*0.5-0.04, label=desc_lab,
             size=2.2, colour=col, hjust=0.5),
    annotate("text", x=(xmin+xmax)/2, y=yb+0.015, label=panel_lab,
             size=2.0, colour="grey40", hjust=0.5)
  )
}
p_eq <- ggplot() + xlim(0,1) + ylim(0,1) + theme_void() +
  theme(plot.background=element_rect(fill="grey98", colour="grey65", linewidth=0.7)) +
  # top: the equation
  annotate("rect", xmin=0.02, xmax=0.98, ymin=0.70, ymax=0.97,
           fill="white", colour="grey75", linewidth=0.5) +
  annotate("text", x=0.5, y=0.84, hjust=0.5, vjust=0.5, size=4.3, parse=TRUE,
           label='gap[g]^{corr}==beta[g]+sum(w[i]*(gap[g*","*i]-beta[g])/sqrt(delta[g*","*i]),i)') +
  annotate("text", x=0.5, y=0.73, hjust=0.5, size=2.6, colour="grey45",
           label=sprintf("per-gene corrected class gap (gene: %s, n=5 training studies)", gene_name)) +
  # term label header
  annotate("text", x=0.5, y=0.66, hjust=0.5, size=2.5, colour="grey35", fontface="italic",
           label="Each term illustrated in a surrounding panel:") +
  # term cards
  mk_card(0.01, 0.30, 0.43, 0.63, tc["gap"],
          'gap[g*","*i]', "per-study class gap", "panels C & F (left)") +
  mk_card(0.32, 0.64, 0.43, 0.63, tc["beta"],
          'beta[g]', "shared pooled estimate", "panel D (above)") +
  mk_card(0.66, 0.99, 0.43, 0.63, tc["delta"],
          'delta[g*","*i]', "within-study scale", "panels E (above) & G (right)") +
  # residual sum card (full width, bottom)
  annotate("rect", xmin=0.01, xmax=0.99, ymin=0.04, ymax=0.39,
           fill=tc["resid"], alpha=0.10, colour=tc["resid"], linewidth=0.4) +
  annotate("text", x=0.5, y=0.27, parse=TRUE, size=3.3, colour=tc["resid"], hjust=0.5,
           label='sum(w[i]*(gap[g*","*i]-beta[g])/sqrt(delta[g*","*i]),i)') +
  annotate("text", x=0.5, y=0.18, size=2.3, colour=tc["resid"], hjust=0.5,
           label="the residual: cancels with equal weights, diverges with real delta") +
  annotate("text", x=0.5, y=0.08, size=2.0, colour="grey40", hjust=0.5,
           label="panels F (left), G (right), H (below)")

# ════════════════════════════════════════════════════════════════════════════
# PANEL G — amplified leftover = (gap - beta) / sqrt(delta)
# ════════════════════════════════════════════════════════════════════════════
d_g <- data.frame(study=sg_fac, val=amplif_gi)
net_sc <- sum(wi_named[as.character(sg_fac)] * amplif_gi)
pG <- ggplot(d_g, aes(val, study, fill=val>0)) +
  geom_vline(xintercept=0, colour="grey55", linewidth=0.35) +
  geom_col(width=0.65, alpha=0.8) +
  annotate("text", x=min(d_g$val)*0.5, y=0.4,
           label=sprintf("weighted sum = %.1f", net_sc),
           hjust=0.5, vjust=0, size=2.5, colour=tc["resid"], fontface="bold") +
  scale_fill_manual(values=c("TRUE"=tc["delta"], "FALSE"="#d6604d"), guide="none") +
  scale_x_continuous(expand=expansion(mult=c(0.12,0.18))) +
  labs(title=sprintf("G  (gap - beta) / sqrt(delta)  [%s]", gene_name),
       subtitle="Same leftovers, now divided by sqrt(delta). Low-scale studies\n(Africa, India, USA) are enormously amplified; they dominate the sum.",
       x="amplified leftover", y=NULL) +
  tag_text("/ sqrt(delta)", tc["delta"]) +
  thl

# ════════════════════════════════════════════════════════════════════════════
# PANEL H — the sum: equal weighting vs real delta (all genes, projected)
# ════════════════════════════════════════════════════════════════════════════
pooled  <- sum(class_effect*e1)
gap_proj <- sapply(seq_len(nb), function(k) sum(gap_i[,k]*e1))
leftover <- gap_proj - pooled
eqw      <- w_i * leftover
scl      <- sapply(seq_len(nb), function(k)
  w_i[k]*sum(interaction[,k]/sqrt(pmax(delta[,k],1e-8))*e1))
levp <- sdsp[order(gap_proj)]; TOT <- "net (sum)"
d_h  <- rbind(
  data.frame(study=sdsp, type="equal\nweighting",        val=eqw),
  data.frame(study=sdsp, type="divide by\nsqrt(scale)", val=scl),
  data.frame(study=TOT,  type="equal\nweighting",        val=sum(eqw)),
  data.frame(study=TOT,  type="divide by\nsqrt(scale)", val=sum(scl)))
d_h$study <- factor(d_h$study, levels=c(TOT, levp))
d_h$type  <- factor(d_h$type,  levels=c("equal\nweighting","divide by\nsqrt(scale)"))
hcol <- c("equal\nweighting"="#9bb8d3","divide by\nsqrt(scale)"=tc["resid"])
pH <- ggplot(d_h, aes(val, study, fill=type)) +
  geom_hline(yintercept=1.5, linetype="dotted", colour="grey60") +
  geom_vline(xintercept=0, colour="grey55", linewidth=0.35) +
  geom_col(position=position_dodge(width=0.7), width=0.62) +
  annotate("text", x=0, y=1.22, hjust=-0.07, size=2.3, colour="#3b5a78",
           label=sprintf("equal: %.1f (cancels)", sum(eqw))) +
  annotate("text", x=sum(scl), y=0.80, hjust=-0.05, size=2.3,
           colour=tc["resid"], fontface="bold",
           label=sprintf("scaled: %.0f (against class)", sum(scl))) +
  scale_fill_manual(values=hcol, name=NULL) +
  scale_x_continuous(expand=expansion(mult=c(0.06,0.18))) +
  labs(title="H  The sum across all genes: equal vs divided by sqrt(scale)",
       subtitle="Projected on the class direction. Equal weighting: per-study contributions cancel.\nDividing by sqrt(scale): low-delta studies dominate -> net residual opposing the class.",
       x="contribution to corrected class gap", y=NULL) +
  th + theme(legend.position="top")

# ════════════════════════════════════════════════════════════════════════════
# PANELS I & J — KNN consequence and fix
# ════════════════════════════════════════════════════════════════════════════
align <- function(ref) {
  comb <- cbind(ref, tst); cb <- c(rep(1L,ncol(ref)), rep(2L,ncol(tst)))
  o <- suppressMessages(ComBat(comb, batch=cb, mod=NULL, ref.batch=1L))
  o[,(ncol(ref)+1):ncol(o), drop=FALSE]
}
hull_pts <- function(d) d[chull(d$x,d$y),]

panelIJ <- function(ref, ttl, sub) {
  tgt <- align(ref); hvg <- order(rowVars(ref), decreasing=TRUE)[1:1000]
  Rh <- ref[hvg,,drop=FALSE]; Th <- tgt[hvg,,drop=FALSE]
  pca <- prcomp(t(Rh), rank.=2); R <- pca$x; Tt <- predict(pca, t(Th))
  dfR <- data.frame(x=R[,1], y=R[,2], cls=clab(lab))
  dfT <- data.frame(x=Tt[,1],y=Tt[,2], cls=clab(ltst))
  hl  <- rbind(cbind(hull_pts(dfR[dfR$cls=="Active TB",]),g="a"),
               cbind(hull_pts(dfR[dfR$cls=="Control",]),  g="c"))
  cen <- rbind(data.frame(x=mean(dfR$x[lab==1]),y=mean(dfR$y[lab==1]),cls="Active TB"),
               data.frame(x=mean(dfR$x[lab==0]),y=mean(dfR$y[lab==0]),cls="Control"))
  set.seed(3); pick <- sample(ncol(Th), 16)
  nn   <- sapply(pick, function(j) which.min(colSums((Rh - Th[,j])^2)))
  links <- data.frame(x=Tt[pick,1],y=Tt[pick,2],xe=R[nn,1],ye=R[nn,2],
                      correct=ifelse(lab[nn]==ltst[pick],"correct","WRONG"))
  pctw <- round(100*mean(links$correct=="WRONG"))
  ggplot() +
    geom_polygon(data=hl, aes(x,y,group=g,fill=ifelse(g=="a","Active TB","Control")), alpha=0.11) +
    geom_point(data=dfR, aes(x,y,colour=cls), size=0.7, alpha=0.28) +
    geom_segment(data=links, aes(x,y,xend=xe,yend=ye,linetype=correct),
                 colour="grey25", linewidth=0.42) +
    geom_point(data=dfT, aes(x,y,colour=cls), shape=17, size=2.1, alpha=0.95) +
    geom_point(data=cen, aes(x,y,fill=cls), shape=23, size=4.0, stroke=1, colour="black") +
    scale_colour_manual(values=ccol, name=NULL, guide=guide_legend(order=1)) +
    scale_fill_manual(values=ccol, guide="none") +
    scale_linetype_manual(values=c("correct"="solid","WRONG"="22"), name=NULL) +
    annotate("text", x=Inf, y=Inf, hjust=1.05, vjust=1.4, size=2.9, fontface="bold",
             colour=ifelse(pctw>50,"#b3261e","#1b7837"),
             label=sprintf("%d%% of shown neighbours: WRONG class", pctw)) +
    labs(title=ttl, subtitle=sub,
         x="PC1 of corrected data (highest variance; what KNN weights most)", y="PC2") +
    th + theme(legend.position="top")
}
pI <- panelIJ(ref_sup,
  "I  Supervised ComBat: held-out samples land by the WRONG class",
  "Triangles = held-out test samples; diamonds = training class centroids.\nActive triangles (red) fall in the training Control cloud (dashed = wrong neighbour).")
pJ <- panelIJ(ref_mo,
  "J  Fix (mean.only): correct class, correct neighbours",
  "Removing only the per-gene scale step (mean.only=TRUE) restores the class axis.\nActive triangles now neighbour training Active samples (solid links).")

# ════════════════════════════════════════════════════════════════════════════
# PANEL K — proportional fix: per-scenario KNN MCC comparison (4 methods)
# ════════════════════════════════════════════════════════════════════════════
prop_fix_fn <- function(dat, bat, lab) {
  ubat <- levels(as.factor(bat)); nb <- length(ubat); N <- ncol(dat)
  ref  <- suppressMessages(ComBat(dat, batch=bat, mod=model.matrix(~lab)))
  des  <- cbind(model.matrix(~ -1+factor(bat)), class=lab)
  B2   <- solve(crossprod(des), crossprod(des, t(dat))); ce2 <- B2["class",]
  nbi2 <- as.numeric(table(bat)); grand2 <- as.numeric(crossprod(nbi2/N, B2[seq_len(nb),]))
  vp2  <- rowMeans((dat - t(des %*% B2))^2)
  sd2  <- (dat - (outer(grand2, rep(1,N)) + outer(ce2, lab))) / sqrt(vp2)
  dh   <- sapply(ubat, function(b) rowVars(sd2[, bat==b, drop=FALSE]))
  out  <- ref
  for (k in seq_along(ubat)) {
    i <- bat == ubat[k]
    factor_gk <- 1/sqrt(pmax(dh[,k], 1e-8)) - 1
    out[,i] <- ref[,i] + outer(ce2 * factor_gk, as.numeric(lab[i]))
  }
  out
}

Rk <- read.csv(file.path(OUT_DIR, "prop_fix_results.csv"))
mlabs <- c(combat="unsupervised\nComBat", combat_sup="supervised\nComBat",
           mean_only="mean.only\n(our fix)", prop_fix="proportional\nscaling fix")
Rk$mlabel <- factor(mlabs[Rk$method], levels=mlabs)
mcol4 <- c("unsupervised\nComBat"="#4e79a7","supervised\nComBat"="#e15759",
           "mean.only\n(our fix)"="#59a14f","proportional\nscaling fix"="#f28e2b")
# medians for annotation
meds <- aggregate(mcc ~ mlabel, Rk, median)

pK <- ggplot(Rk, aes(mlabel, mcc, colour=mlabel)) +
  geom_hline(yintercept=0, colour="grey50", linewidth=0.4) +
  geom_jitter(width=0.18, size=1.5, alpha=0.7) +
  stat_summary(fun=median, geom="crossbar", width=0.55, fatten=2.5, colour="black", linewidth=0.6) +
  scale_colour_manual(values=mcol4, guide="none") +
  labs(title="K  Proposed proportional scaling fix: empirical KNN MCC (24 scenarios)",
       subtitle="Each dot = one leave-one-study-out scenario (n=2..5 training studies). Bar = median.\nProp fix prevents most sign flips (2/24 below chance) but substantially underperforms mean.only.",
       x=NULL, y="KNN MCC (held-out)") +
  th

# ════════════════════════════════════════════════════════════════════════════
# PANEL L — what each method does to gene RPS4Y1's corrected class gap
# ════════════════════════════════════════════════════════════════════════════
ref_prop <- prop_fix_fn(dat, bat, lab)
wvec     <- function(m, l) rowMeans(m[,l==1,drop=FALSE]) - rowMeans(m[,l==0,drop=FALSE])
ref_unsup <- suppressMessages(ComBat(dat, batch=bat, mod=NULL))

gaps <- data.frame(
  method = factor(c("raw\n(before ComBat)",
                    "unsupervised\nComBat",
                    "mean.only\n(our fix)",
                    "proportional\nscaling fix",
                    "supervised\nComBat"),
                  levels=c("raw\n(before ComBat)","unsupervised\nComBat",
                           "mean.only\n(our fix)","proportional\nscaling fix",
                           "supervised\nComBat")),
  gap = c(wvec(dat,     lab)[gi],
          wvec(ref_unsup,lab)[gi],
          wvec(ref_mo,   lab)[gi],
          wvec(ref_prop, lab)[gi],
          wvec(ref_sup,  lab)[gi])
)
gcol5 <- c("raw\n(before ComBat)"="grey50",
           "unsupervised\nComBat"="#4e79a7",
           "mean.only\n(our fix)"="#59a14f",
           "proportional\nscaling fix"="#f28e2b",
           "supervised\nComBat"="#e15759")

pL <- ggplot(gaps, aes(gap, method, fill=method)) +
  geom_vline(xintercept=0, colour="grey55", linewidth=0.35) +
  geom_col(width=0.65) +
  geom_text(aes(label=sprintf("%+.3f", gap),
                hjust=ifelse(gap>=0, -0.1, 1.1)),
            size=2.8, colour="grey20") +
  scale_fill_manual(values=gcol5, guide="none") +
  scale_x_continuous(expand=expansion(mult=c(0.18,0.18))) +
  labs(title=sprintf("L  Corrected class gap for gene %s under each method", gene_name),
       subtitle=paste0("Active minus Control mean expression after correction. Positive = Active > Control (correct).\n",
                       "Supervised ComBat flips the sign; both fixes restore it, but prop fix overamplifies the gap."),
       x=sprintf("corrected Active-Control gap  [%s, log2]", gene_name), y=NULL) +
  th

# ════════════════════════════════════════════════════════════════════════════
# ASSEMBLE
# ════════════════════════════════════════════════════════════════════════════
layout_matrix <- rbind(
  c(1,1,1, 2,2,2),     # row 1: A | B
  c(3,3, 4,4, 5,5),    # row 2: C | D | E
  c(6,6, 7,7, 8,8),    # row 3: F | EQ | G
  c(9,9, 10,10, 11,11),# row 4: H | I  | J
  c(12,12,12, 13,13,13) # row 5: K | L
)
panels <- list(pA, pB, pC, pD, pE, pF, p_eq, pG, pH, pI, pJ, pK, pL)
g <- arrangeGrob(grobs=panels, layout_matrix=layout_matrix,
                 heights=c(1.0, 1.4, 1.6, 1.6, 1.3))
ggsave(file.path(OUT_DIR,"equation_figure.png"), g, width=17, height=27, dpi=180)
ggsave(file.path(OUT_DIR,"equation_figure.pdf"), g, width=17, height=27)
cat("-> equation_figure.png/.pdf\n")
