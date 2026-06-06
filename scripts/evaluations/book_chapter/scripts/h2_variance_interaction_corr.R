#!/usr/bin/env Rscript
# h2_variance_interaction_corr.R — Does option2's within-class variance decouple
# from the class-by-batch interaction, while ComBat's delta picks it up?
# H2 prediction: cor(vl_b, |gap_b - beta|) ≈ 0; cor(delta.star_b, |gap_b - beta|) > 0.

suppressMessages({ library(sva); library(matrixStats) })
setwd("/home/phr23/confounded_analysis/scripts/evaluations/book_chapter")
load("data/TB_real_data.RData")
OUT_DIR <- "outputs/diagnostics/hypothesis_tests"

ALL_STUDIES <- c("GSE37250_SA","USA","India","GSE37250_M","Africa","GSE39941_M")
bin <- function(l) as.integer(ifelse(l %in% c("1",1,"Active"), 1L, 0L))
train_studies <- function(n, test) ALL_STUDIES[ALL_STUDIES!=test][seq_len(n)]
log_safe <- function(m) { if(max(m,na.rm=TRUE)>100){ mn<-min(m,na.rm=TRUE); if(mn<0) m<-m-mn; m<-log2(m+1) }; m }

aprior <- function(gg) { m<-mean(gg); s2<-var(gg); (2*s2+m^2)/s2 }
bprior <- function(gg) { m<-mean(gg); s2<-var(gg); (m*s2+m^3)/s2 }
postvar_fn <- function(sum2, n, a, b) (0.5*sum2+b)/(n/2+a-1)

rows <- list()
for(nn in 2:5) for(tst in ALL_STUDIES) {
  refs <- train_studies(nn, tst)
  cg <- Reduce(intersect, lapply(c(refs,tst), function(s) rownames(dat_lst[[s]])))
  dat <- log_safe(do.call(cbind, lapply(refs, function(s) dat_lst[[s]][cg,])))
  bat <- as.factor(do.call(c, lapply(refs, function(s) rep(s, ncol(dat_lst[[s]])))))
  lab <- do.call(c, lapply(refs, function(s) bin(label_lst[[s]])))
  ubat <- levels(bat); N <- ncol(dat)

  # Global LS fit (supervised)
  des <- cbind(model.matrix(~-1+bat), class=lab)
  B <- solve(crossprod(des), crossprod(des, t(dat)))
  ce <- B["class",]
  nbi <- as.numeric(table(bat)); grand <- as.numeric(crossprod(nbi/N, B[seq_len(length(ubat)),]))
  vp_ls <- rowMeans((dat - t(des %*% B))^2)
  sm <- outer(grand, rep(1,N)) + outer(ce, lab)
  Z <- (dat - sm) / sqrt(pmax(vp_ls, 1e-8))

  for(k in seq_along(ubat)) {
    b <- ubat[k]; i <- bat==b; n_b <- sum(i)
    n1 <- sum(i & lab==1); n0 <- sum(i & lab==0)
    if(n1 < 1 || n0 < 1) next

    # Per-study raw class gap
    gap_b <- rowMeans(dat[,i & lab==1,drop=FALSE]) - rowMeans(dat[,i & lab==0,drop=FALSE])
    # Interaction = gap_b minus pooled class effect
    interaction <- abs(gap_b - ce)

    # Option 2: within-class pooled variance
    v1 <- if(n1>1) rowVars(dat[,i & lab==1,drop=FALSE]) else rep(0, nrow(dat))
    v0 <- if(n0>1) rowVars(dat[,i & lab==0,drop=FALSE]) else rep(0, nrow(dat))
    vl_b <- (pmax(n1-1,1)*v1 + pmax(n0-1,1)*v0) / pmax(n1+n0-2, 1)

    # ComBat delta_raw (same as option3 dhat)
    Zi <- Z[,i,drop=FALSE]; ghat_b <- rowMeans(Zi)
    dr_b <- rowVars(Zi - ghat_b)

    # ComBat delta.star (EB-shrunk)
    if(n_b > 2 && var(dr_b) > 1e-10) {
      a_b <- aprior(dr_b); s_b <- bprior(dr_b)
      ds_b <- postvar_fn(dr_b*(n_b-1), n_b, a_b, s_b)
    } else {
      ds_b <- dr_b
    }

    # Summarize: median across genes (representative of bulk behavior)
    rows[[length(rows)+1]] <- data.frame(
      n=nn, test=tst, batch=b,
      # Spearman correlations: interaction vs variance estimate
      cor_vl_interaction    = cor(vl_b, interaction, method="spearman", use="complete.obs"),
      cor_dr_interaction    = cor(dr_b, interaction, method="spearman", use="complete.obs"),
      cor_dstar_interaction = cor(ds_b, interaction, method="spearman", use="complete.obs"),
      # Medians for context
      med_vl   = median(vl_b),
      med_dr   = median(dr_b),
      med_dstar= median(ds_b),
      med_interaction = median(interaction),
      n_b=n_b, stringsAsFactors=FALSE)
  }
  cat("done n=",nn," test=",tst,"\n",sep="")
}
res <- do.call(rbind, rows)
write.csv(res, file.path(OUT_DIR, "h2_variance_interaction_corr.csv"), row.names=FALSE)

cat("\n=== H2 Summary: Spearman cor(variance_estimate, |gap - beta|) ===\n")
cat("Method             median cor  mean cor\n")
cat("Option2 (vl)      ", round(median(res$cor_vl_interaction),3),
    "       ", round(mean(res$cor_vl_interaction),3), "\n")
cat("delta_raw (opt3)  ", round(median(res$cor_dr_interaction),3),
    "       ", round(mean(res$cor_dr_interaction),3), "\n")
cat("delta.star (ComBat)", round(median(res$cor_dstar_interaction),3),
    "       ", round(mean(res$cor_dstar_interaction),3), "\n")
cat("\nH2 verdict: if vl correlation is near 0 but delta.star > 0, H2 is confirmed.\n")
cat("-> h2_variance_interaction_corr.csv\n")
