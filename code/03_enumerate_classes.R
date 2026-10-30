# RQ1: class enumeration for the four AI-use purpose items (ST438Q01-04).
# Weighted latent class models with 1-8 classes (500 random starts each), bivariate
# residuals, a check of the estimator against poLCA, and replication in two random
# splits of schools into halves.

source("code/common.R")
library(poLCA)

set.seed(20261029)
J <- length(items); K_max <- 8; n_start <- 500; tol <- 1e-11; maxit <- 20000

d <- load_sample()
n <- nrow(d)
y <- as.matrix(d[, items])
w <- d$W_FSTUWT / mean(d$W_FSTUWT)            # rescaled to sum to n
cat("students", n, " schools", length(unique(d$CNTSCHID)), "\n")

pat <- collapse_patterns(y, w)
oh <- one_hot(pat$y, R)
cat("unique response patterns", length(pat$w), "\n")

# ---- Enumeration (weighted) ----
fits <- list(); stats <- list()
for (Kc in 1:K_max) {
  t0 <- Sys.time()
  f <- lca_fit_multi(oh, pat$w, Kc, R, nstart = if (Kc == 1) 1 else n_start,
                     tol = tol, maxit = maxit)
  f <- lca_order(f, R)
  fits[[Kc]] <- f
  stats[[Kc]] <- lca_fit_stats(f, pat$w, n, J, R)
  cat(sprintf("K=%d  LL=%.2f  best LL replicated %d/%d  (%.0fs)\n", Kc, f$ll,
              f$best_replicated, if (Kc == 1) 1 else n_start,
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}
fit_tab <- do.call(rbind, stats)

# Bivariate residuals: Pearson chi-square of each item pair divided by its 16 df
bvr_pairs <- function(f) {
  b <- c()
  for (j in 1:(J - 1)) for (l in (j + 1):J) {
    ok <- !is.na(pat$y[, j]) & !is.na(pat$y[, l])
    obs <- xtabs(pat$w[ok] ~ factor(pat$y[ok, j], 1:R) + factor(pat$y[ok, l], 1:R))
    exp_n <- sum(obs) * Reduce(`+`, lapply(seq_len(f$K), function(k)
      f$pi[k] * outer(f$rho[[j]][k, ], f$rho[[l]][k, ])))
    b[paste(items[j], items[l], sep = " x ")] <- sum((obs - exp_n)^2 / exp_n) / (R - 1)^2
  }
  b
}
fit_tab$max_BVR <- sapply(fits, function(f) max(bvr_pairs(f)))
print(format(fit_tab, digits = 4), row.names = FALSE)
cat("\nBivariate residuals, six-class model:\n"); print(round(bvr_pairs(fits[[K]]), 1))
write.csv(fit_tab, res_path("lca_fit_indices.csv"), row.names = FALSE)

prof <- do.call(rbind, lapply(fits[-1], lca_profile_table, items = items))
write.csv(prof, res_path("lca_item_probabilities.csv"), row.names = FALSE)

# Within-person shape of the raw answers (complete cases, weighted shares)
cc <- complete.cases(y); yc <- y[cc, ]; wc <- w[cc]
wshare <- function(x) round(weighted.mean(x, wc), 3)
rng <- apply(yc, 1, function(r) max(r) - min(r))
cat("\nRaw answers (weighted, complete cases):\n")
cat("  all never", wshare(rowSums(yc == 1) == 4), " all daily", wshare(rowSums(yc == 5) == 4),
    " identical other", wshare(rng == 0 & yc[, 1] %in% 2:4), "\n")
cat("  summarize or draft above learn by >= 1 level", wshare(pmax(yc[, 1], yc[, 3]) - yc[, 4] >= 1),
    "; by >= 2 levels", wshare(pmax(yc[, 1], yc[, 3]) - yc[, 4] >= 2), "\n")
cat("  within-person range 0..4:", round(tapply(wc, rng, sum) / sum(wc), 3), "\n")

cat("\nClass shares and expected response level (1 = never ... 5 = daily):\n")
for (Kc in 3:K_max) {
  f <- fits[[Kc]]
  lev <- sapply(f$rho, function(r) r %*% seq_len(R))
  tab <- data.frame(class = seq_len(Kc), share = round(f$pi, 3), round(lev, 2))
  names(tab)[3:6] <- c("summarize", "research", "draft", "learn")
  cat("K =", Kc, "\n"); print(tab, row.names = FALSE)
}

# ---- Estimator check: unweighted EM against poLCA ----
pat_u <- collapse_patterns(y, rep(1, n))
oh_u <- one_hot(pat_u$y, R)
pd <- as.data.frame(y); names(pd) <- paste0("u", 1:J)
check <- data.frame()
for (Kc in 2:5) {
  own <- lca_fit_multi(oh_u, pat_u$w, Kc, R, nstart = 30)
  pl <- poLCA(cbind(u1, u2, u3, u4) ~ 1, data = pd, nclass = Kc, nrep = 30,
              na.rm = FALSE, verbose = FALSE, maxiter = 5000, tol = 1e-10)
  check <- rbind(check, data.frame(K = Kc, LL_own = own$ll, LL_poLCA = pl$llik,
                                   diff = own$ll - pl$llik))
}
cat("\nUnweighted log-likelihood, own EM vs poLCA:\n"); print(check, row.names = FALSE)
write.csv(check, res_path("lca_check_polca.csv"), row.names = FALSE)

# ---- Replication in two random splits of schools into halves ----
schools <- unique(d$CNTSCHID)
split <- data.frame(); halves <- list()
for (s in 1:2) {
  half_a <- sample(schools, floor(length(schools) / 2))
  halves[[s]] <- half_a
  for (h in c("A", "B")) {
    idx <- if (h == "A") d$CNTSCHID %in% half_a else !(d$CNTSCHID %in% half_a)
    wh <- d$W_FSTUWT[idx] / mean(d$W_FSTUWT[idx])
    ph <- collapse_patterns(y[idx, ], wh)
    ohh <- one_hot(ph$y, R)
    for (Kc in 2:K_max) {
      fh <- lca_fit_multi(ohh, ph$w, Kc, R, nstart = 100, tol = tol, maxit = maxit)
      m <- lca_match(fits[[Kc]], fh)
      st <- lca_fit_stats(fh, ph$w, sum(idx), J, R)
      split <- rbind(split, data.frame(split = s, half = h, n = sum(idx), K = Kc, BIC = st$BIC,
                                       aBIC = st$aBIC, min_class = st$min_class,
                                       mad_vs_full = m$mad,
                                       shares = paste(round(m$fit$pi, 3), collapse = "/")))
    }
  }
}
cat("\nSplit-half replication:\n"); print(split, row.names = FALSE)
write.csv(split, res_path("lca_split_half.csv"), row.names = FALSE)

saveRDS(list(items = items, fits = lapply(fits, function(f) f[c("K", "pi", "rho", "ll")]),
             fit_tab = fit_tab, half_a_schools = halves),
        res_path("lca_fits.rds"))
