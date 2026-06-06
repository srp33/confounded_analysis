#!/usr/bin/env Rscript
# h1_delta_comparison.R — Compare option3 delta_raw vs ComBat EB-shrunk delta.star
# H1: option3 fails because raw (unshrunken) dhat contains very small values that
#     standard ComBat's EB shrinkage would have pulled up, causing more extreme amplification.

suppressMessages({ library(sva); library(matrixStats); library(ggplot2); library(gridExtra) })
setwd("/home/phr23/confounded_analysis/scripts/evaluations/book_chapter")
load("data/TB_real_data.RData")
OUT_DIR <- "outputs/diagnostics/hypothesis_tests"

ALL_STUDIES <- c("GSE37250_SA","USA","India","GSE37250_M","Africa","GSE39941_M")
bin <- function(l) as.integer(ifelse(l %in% c("1",1,"Active"), 1L, 0L))
train_studies <- function(n, test) ALL_STUDIES[ALL_STUDIES!=test][seq_len(n)]
log_safe <- function(m) { if(max(m,na.rm=TRUE)>100){ mn<-min(m,na.rm=TRUE); if(mn<0) m<-m-mn; m<-log2(m+1) }; m }

aprior <- function(gg) { m <- mean(gg); s2 <- var(gg); (2*s2 + m^2)/s2 }
bprior <- function(gg) { m <- mean(gg); s2 <- var(gg); (m*s2 + m^3)/s2 }
postvar_fn <- function(sum2, n, a, b) (0.5*sum2 + b) / (n/2 + a - 1)

compute_deltas <- function(dat, bat, lab) {
  ubat <- levels(bat); nb_arr <- as.numeric(table(bat)); N <- ncol(dat)
  des <- cbind(model.matrix(~-1+bat), class=lab)
  B <- solve(crossprod(des), crossprod(des, t(dat)))
  ce <- B["class",]; nbi <- nb_arr
  grand <- as.numeric(crossprod(nbi/N, B[seq_len(length(ubat)),]))
  vp <- rowMeans((dat - t(des %*% B))^2)
  sm <- outer(grand, rep(1,N)) + outer(ce, lab)
  Z <- (dat - sm) / sqrt(pmax(vp, 1e-8))

  delta_raw  <- matrix(NA, nrow(dat), length(ubat), dimnames=list(NULL, ubat))
  delta_star <- matrix(NA, nrow(dat), length(ubat), dimnames=list(NULL, ubat))
  for(k in seq_along(ubat)) {
    i <- bat == ubat[k]; n_b <- sum(i)
    Zi <- Z[,i,drop=FALSE]
    ghat <- rowMeans(Zi)
    dr <- rowVars(Zi - ghat)
    delta_raw[,k] <- dr
    if(n_b > 2 && var(dr) > 1e-10) {
      a_b <- aprior(dr); s_b <- bprior(dr)
      sum_sq <- dr * (n_b - 1)
      ds <- postvar_fn(sum_sq, n_b, a_b, s_b)
      delta_star[,k] <- pmax(ds, 1e-8)
    } else {
      delta_star[,k] <- dr
    }
  }
  list(delta_raw=delta_raw, delta_star=delta_star, beta=ce, vp=vp, Z=Z, sm=sm, bat=bat, lab=lab, ubat=ubat)
}

# ── Focal scenario: n=5, test=GSE37250_SA ──────────────────────────────────────
n <- 5; test <- "GSE37250_SA"; refs <- train_studies(n, test)
cg <- Reduce(intersect, lapply(c(refs,test), function(s) rownames(dat_lst[[s]])))
dat <- log_safe(do.call(cbind, lapply(refs, function(s) dat_lst[[s]][cg,])))
bat <- as.factor(do.call(c, lapply(refs, function(s) rep(s, ncol(dat_lst[[s]])))))
lab <- do.call(c, lapply(refs, function(s) bin(label_lst[[s]])))
out <- compute_deltas(dat, bat, lab)

# ── Per-gene, per-batch: flip thresholds ───────────────────────────────────────
# For gene g, batch b: flip when delta < flip_thr = max(0, 1 - gap_raw/beta)
# For ComBat (delta.star): flip_thr; for option3 (delta_raw): flip_thr^2
gaps <- sapply(out$ubat, function(b) {
  i <- bat == b
  rowMeans(dat[,i & lab==1,drop=FALSE]) - rowMeans(dat[,i & lab==0,drop=FALSE])
})  # G x B matrix of raw per-study class gaps
beta <- out$beta

# Flip threshold per gene per batch (only defined when gap_raw < beta and beta > 0)
flip_thr <- matrix(NA, nrow(dat), length(out$ubat))
for(k in seq_along(out$ubat)) {
  b_g <- beta; g_g <- gaps[,k]
  # Only valid when beta > 0 and gap_raw < beta
  valid <- b_g > 0 & g_g < b_g
  flip_thr[valid, k] <- pmax(0, 1 - g_g[valid]/b_g[valid])
}

# Count flipped genes per batch under each method
flip_combat <- colSums(out$delta_star < flip_thr, na.rm=TRUE)
flip_opt3   <- colSums(out$delta_raw  < flip_thr^2, na.rm=TRUE)  # option3 needs delta < thr^2

cat("=== Focal scenario n=5, test=GSE37250_SA ===\n")
cat("Per-batch flip counts (genes with inverted class axis):\n")
ft <- data.frame(batch=out$ubat, flip_combat=flip_combat, flip_opt3=flip_opt3,
                 n_b=as.numeric(table(bat)), stringsAsFactors=FALSE)
print(ft); cat("\n")

# Summary stats on delta distributions
dr_vec  <- as.vector(out$delta_raw[, , drop=FALSE])
ds_vec  <- as.vector(out$delta_star[, , drop=FALSE])
cat("delta_raw: median=", round(median(dr_vec,na.rm=TRUE),3),
    " pct<1=", round(mean(dr_vec<1,na.rm=TRUE)*100,1), "%",
    " pct<0.1=", round(mean(dr_vec<0.1,na.rm=TRUE)*100,1), "%\n")
cat("delta_star: median=", round(median(ds_vec,na.rm=TRUE),3),
    " pct<1=", round(mean(ds_vec<1,na.rm=TRUE)*100,1), "%",
    " pct<0.1=", round(mean(ds_vec<0.1,na.rm=TRUE)*100,1), "%\n\n")

# ── Across all 24 scenarios: fraction of flipped genes ────────────────────────
rows <- list()
for(nn in 2:5) for(tst in ALL_STUDIES) {
  refs2 <- train_studies(nn, tst)
  cg2 <- Reduce(intersect, lapply(c(refs2,tst), function(s) rownames(dat_lst[[s]])))
  dat2 <- log_safe(do.call(cbind, lapply(refs2, function(s) dat_lst[[s]][cg2,])))
  bat2 <- as.factor(do.call(c, lapply(refs2, function(s) rep(s, ncol(dat_lst[[s]])))))
  lab2 <- do.call(c, lapply(refs2, function(s) bin(label_lst[[s]])))
  o2 <- compute_deltas(dat2, bat2, lab2)
  gaps2 <- sapply(o2$ubat, function(b) {
    i <- bat2==b; rowMeans(dat2[,i&lab2==1,drop=FALSE]) - rowMeans(dat2[,i&lab2==0,drop=FALSE])
  })
  ft2 <- matrix(NA, nrow(dat2), length(o2$ubat))
  for(k in seq_along(o2$ubat)) {
    bg <- o2$beta; gg <- gaps2[,k]; valid <- bg>0 & gg<bg
    ft2[valid,k] <- pmax(0, 1 - gg[valid]/bg[valid])
  }
  fc <- sum(o2$delta_star < ft2, na.rm=TRUE)
  fo <- sum(o2$delta_raw  < ft2^2, na.rm=TRUE)
  total_valid <- sum(!is.na(ft2))
  rows[[length(rows)+1]] <- data.frame(n=nn, test=tst,
    flip_combat_n=fc, flip_opt3_n=fo, total_valid=total_valid,
    flip_combat_pct=100*fc/max(total_valid,1), flip_opt3_pct=100*fo/max(total_valid,1),
    median_dstar=median(o2$delta_star,na.rm=TRUE),
    median_draw=median(o2$delta_raw,na.rm=TRUE),
    stringsAsFactors=FALSE)
  cat("done n=",nn," test=",tst,"\n",sep="")
}
res <- do.call(rbind, rows)
write.csv(res, file.path(OUT_DIR, "h1_delta_comparison.csv"), row.names=FALSE)

cat("\n=== Summary across 24 scenarios ===\n")
cat("Median delta_raw:", round(median(res$median_draw),3), "\n")
cat("Median delta_star:", round(median(res$median_dstar),3), "\n")
cat("Ratio flip_opt3/flip_combat: min=",round(min(res$flip_opt3_pct/pmax(res$flip_combat_pct,0.01)),2),
    " median=",round(median(res$flip_opt3_pct/pmax(res$flip_combat_pct,0.01)),2),
    " max=",round(max(res$flip_opt3_pct/pmax(res$flip_combat_pct,0.01)),2),"\n")

# ── Plot: delta_raw vs delta.star for focal scenario ──────────────────────────
dr_all <- data.frame(value=dr_vec, type="delta_raw (option3)")
ds_all <- data.frame(value=ds_vec, type="delta.star (ComBat EB)")
dens_dat <- rbind(dr_all, ds_all)
dens_dat <- dens_dat[dens_dat$value < 10,]  # trim extreme tail for display
p1 <- ggplot(dens_dat, aes(value, fill=type, colour=type)) +
  geom_density(alpha=0.35) + scale_x_log10() +
  geom_vline(xintercept=1, linetype="dashed", colour="grey40") +
  labs(title="H1: delta_raw (option3) vs delta.star (ComBat EB)",
       subtitle="Focal n=5 test=GSE37250_SA. Values <1 amplify the residual; small values cause flips.",
       x="delta value (log scale)", y="density") +
  theme_minimal(base_size=10)

# Scatter: delta_raw vs delta.star (per-gene per-batch)
scat <- data.frame(delta_raw=as.vector(out$delta_raw), delta_star=as.vector(out$delta_star))
scat <- scat[complete.cases(scat) & scat$delta_raw<10 & scat$delta_star<10,]
p2 <- ggplot(scat[sample(nrow(scat), min(5000,nrow(scat))),], aes(delta_raw, delta_star)) +
  geom_point(alpha=0.15, size=0.4) + geom_abline(slope=1, intercept=0, colour="red") +
  geom_hline(yintercept=1, linetype="dashed", colour="grey50") +
  geom_vline(xintercept=1, linetype="dashed", colour="grey50") +
  scale_x_log10() + scale_y_log10() +
  labs(title="delta_raw vs delta.star per gene per batch",
       subtitle="Points below y=1 (dashed): ComBat risk zone; below red line: EB shrank delta down",
       x="delta_raw (option3)", y="delta.star (ComBat EB)") +
  theme_minimal(base_size=10)

g <- arrangeGrob(p1, p2, ncol=2)
ggsave(file.path(OUT_DIR,"h1_delta_comparison.png"), g, width=12, height=5, dpi=150)
cat("-> h1_delta_comparison.csv + .png\n")
