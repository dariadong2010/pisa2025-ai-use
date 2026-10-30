# RQ4: class differences in response style by item set, and the share of the
# daily-vs-non-user disposition difference that remains after adjusting for each
# version of the ERS index.
# Each disposition index is standardised once with its weighted SD; adjusted versions are
# the residuals of that z-score on an ERS index and balanced acquiescence. Raw and adjusted
# composites are divided by the weighted SD of the raw composite (1 = one SD of the
# observed composite). SEs: LCA, D and BCH repeated for the 80 BRR-Fay replicates.

source("code/common.R")

ref <- load_reference()
d <- load_sample(response_style = TRUE)
y <- as.matrix(d[, items])
oh_all <- one_hot(y, R)
w0 <- d$W_FSTUWT

wmean <- function(x, w) { ok <- !is.na(x); sum(w[ok] * x[ok]) / sum(w[ok]) }
wsd <- function(x, w) { ok <- !is.na(x); m <- wmean(x, w); sqrt(sum(w[ok] * (x[ok] - m)^2) / sum(w[ok])) }
disp <- c("CURIO", "PERSEV", "COGABIL", "ENPROBS")
Z <- sapply(disp, function(v) (d[[v]] - wmean(d[[v]], w0)) / wsd(d[[v]], w0))
comp_raw <- ifelse(rowSums(!is.na(Z)) >= 3, rowMeans(Z, na.rm = TRUE), NA)
sd_comp <- wsd(comp_raw, w0)
cat("weighted SD of the raw composite:", round(sd_comp, 3), "\n")
d$DISP <- comp_raw / sd_comp

variants <- c("ERS", "ERS_pos", "ERS_neg", "ERS_neu", "ERS_norel", "ERS_ict")
cat("\nR2 of each disposition index on the ERS version and balanced acquiescence:\n")
for (v in variants) {
  Rz <- sapply(seq_along(disp), function(j) {
    f <- lm(Z[, j] ~ d[[v]] + d$ARS_bal, weights = w0, na.action = na.exclude)
    if (j == 1) cat(sprintf("  %-10s", v))
    cat(sprintf(" %s %.3f", disp[j], summary(f)$r.squared))
    residuals(f)
  })
  cat("\n")
  d[[paste0("DISP_adj_", v)]] <- ifelse(rowSums(!is.na(Rz)) >= 3, rowMeans(Rz, na.rm = TRUE), NA) / sd_comp
}

outcomes <- c(variants, "EXT_agree", "EXT_disagree", "ARS_bal", "DISP", paste0("DISP_adj_", variants))

run <- function(wraw) {
  w <- wraw / mean(wraw)
  pat <- collapse_patterns(y, w)
  f <- lca_em(one_hot(pat$y, R), pat$w, K, R, start = ref, tol = 1e-11, maxit = 20000)
  f <- lca_match(ref, f)$fit
  post <- lca_posterior(f, oh_all)
  modal <- max.col(post, ties.method = "first")
  D <- t(sapply(seq_len(K), function(k)
    sapply(seq_len(K), function(s) sum(w * post[, k] * (modal == s)) / sum(w * post[, k]))))
  H <- solve(D)[modal, ]
  sapply(outcomes, function(o) {
    z <- d[[o]]; ok <- !is.na(z)
    colSums(w[ok] * H[ok, ] * z[ok]) / colSums(w[ok] * H[ok, ])
  })                                                    # K x outcomes
}

est <- run(w0)
reps <- lapply(1:80, function(r) {
  if (r %% 20 == 0) cat("replicate", r, "\n")
  run(d[[paste0("W_FSTURWT", r)]])
})
se_of <- function(g) {
  q <- g(est)
  qr <- matrix(sapply(reps, g), ncol = 80)
  sqrt(rowSums((qr - as.vector(q))^2) / (80 * 0.25))
}

mean_tab <- do.call(rbind, lapply(outcomes, function(o)
  data.frame(outcome = o, class = class_names, mean = est[, o], se = se_of(function(e) e[, o]))))
write.csv(mean_tab, res_path("rs_class_means.csv"), row.names = FALSE)

# Daily all-purpose users against every other class
con_tab <- do.call(rbind, lapply(outcomes, function(o) {
  do.call(rbind, lapply(1:5, function(b) {
    g <- function(e) e[6, o] - e[b, o]
    dd <- g(est); s <- se_of(g)
    data.frame(outcome = o, contrast = paste("daily all-purpose vs", class_names[b]),
               diff = dd, se = s, p = 2 * pnorm(-abs(dd / s)))
  }))
}))
con_tab$p_holm <- ave(con_tab$p, con_tab$outcome, FUN = function(p) p.adjust(p, "holm"))
write.csv(con_tab, res_path("rs_contrasts.csv"), row.names = FALSE)

# Share of the daily-vs-non-user disposition difference left after each adjustment
rem_tab <- do.call(rbind, lapply(variants, function(v) {
  a <- paste0("DISP_adj_", v)
  g_raw <- function(e) e[6, "DISP"] - e[1, "DISP"]
  g_adj <- function(e) e[6, a] - e[1, a]
  g_ratio <- function(e) g_adj(e) / g_raw(e)
  data.frame(ers_version = v, diff_raw = g_raw(est), se_raw = se_of(g_raw),
             diff_adj = g_adj(est), se_adj = se_of(g_adj),
             share_remaining = g_ratio(est), se_share = se_of(g_ratio))
}))
write.csv(rem_tab, res_path("rs_remaining.csv"), row.names = FALSE)

cat("\nClass means (SE):\n")
m <- matrix(sprintf("%.3f (%.3f)", mean_tab$mean, mean_tab$se), nrow = K,
            dimnames = list(class_names, outcomes))
print(noquote(t(m)))
cat("\nDaily all-purpose vs non-users and vs frequent all-purpose:\n")
show <- con_tab[con_tab$contrast %in% c("daily all-purpose vs non-user", "daily all-purpose vs frequent all-purpose"), ]
show[, c("diff", "se")] <- round(show[, c("diff", "se")], 3); show$p_holm <- signif(show$p_holm, 2)
print(show[, c("outcome", "contrast", "diff", "se", "p_holm")], row.names = FALSE)
cat("\nDisposition difference (daily - non-user) remaining after adjustment:\n")
print(format(rem_tab, digits = 3), row.names = FALSE)
