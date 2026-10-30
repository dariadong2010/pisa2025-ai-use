# Robustness of the class solution:
#   A. unweighted six-class LCA compared with the weighted solution;
#   B. unweighted BCH means for dispositions, science and their difference (point estimates);
#   C. multilevel LCA (multilevLCA, unweighted; students within schools) with 2-4
#      school-level classes, compared with the single-level solution. The package's
#      default start values are used.

source("code/common.R")
library(multilevLCA)

set.seed(20261102)
ref <- load_reference()
levels_of <- function(f) round(sapply(f$rho, function(r) r %*% seq_len(R)), 2)

d <- load_sample()
y <- as.matrix(d[, items])
n <- nrow(d)

# ---- A. Unweighted LCA ----
pat_u <- collapse_patterns(y, rep(1, n))
fu <- lca_fit_multi(one_hot(pat_u$y, R), pat_u$w, K, R, nstart = 200, tol = 1e-11, maxit = 20000)
mu <- lca_match(ref, fu)
fu <- mu$fit
cat("A. Unweighted six-class LCA: best LL replicated", fu$best_replicated, "/ 200\n")
cat("   mean abs. difference in item probabilities vs weighted:", round(mu$mad, 4), "\n")
cat("   unweighted BIC:", round(-2 * fu$ll + (K - 1 + K * length(items) * (R - 1)) * log(n), 1), "\n")
print(data.frame(class = class_names, weighted = round(ref$pi, 3), unweighted = round(fu$pi, 3),
                 levels_of(fu)), row.names = FALSE)

# ---- B. Unweighted BCH (point estimates) ----
bch_H <- function(f, w) {
  post <- lca_posterior(f, one_hot(y, R))
  modal <- max.col(post, ties.method = "first")
  D <- t(sapply(seq_len(K), function(k)
    sapply(seq_len(K), function(s) sum(w * post[, k] * (modal == s)) / sum(w * post[, k]))))
  solve(D)[modal, ]
}
bch_mean <- function(z, w, H) {
  ok <- !is.na(z)
  colSums(w[ok] * H[ok, ] * z[ok]) / colSums(w[ok] * H[ok, ])
}
wz <- d$W_FSTUWT
zstd <- function(x) {
  ok <- !is.na(x); m <- weighted.mean(x[ok], wz[ok])
  (x - m) / sqrt(weighted.mean((x[ok] - m)^2, wz[ok]))
}
zd <- sapply(c("CURIO", "PERSEV", "COGABIL", "ENPROBS"), function(v) zstd(d[[v]]))
disp_z <- ifelse(rowSums(!is.na(zd)) >= 3, rowMeans(zd, na.rm = TRUE), NA)
outcome_means <- function(f, w) {
  H <- bch_H(f, w)
  sci <- rowMeans(sapply(1:M, function(m) bch_mean(d[[paste0("PV", m, "SCIE")]], w, H)))
  gap <- rowMeans(sapply(1:M, function(m) bch_mean(disp_z - zstd(d[[paste0("PV", m, "SCIE")]]), w, H)))
  data.frame(class = class_names, DISP_z = round(bch_mean(disp_z, w, H), 2),
             SCIE = round(sci), GAP_SCIE = round(gap, 2))
}
bw <- outcome_means(ref, wz / mean(wz))
bu <- outcome_means(fu, rep(1, n))
cat("\nB. BCH means, weighted (left) vs unweighted (right):\n")
print(cbind(bw, bu[, -1]), row.names = FALSE)

# ---- C. Multilevel LCA (students in schools) ----
dd <- data.frame(y - 1, school = d$CNTSCHID)
names(dd)[1:4] <- items
to_fit <- function(m) {
  rho <- lapply(seq_along(items), function(j) t(m$mPhi[(j - 1) * R + seq_len(R), , drop = FALSE]))
  pi_low <- if (is.null(m$vOmega)) as.vector(m$vPi) else as.vector(m$mPi %*% m$vOmega)
  list(K = K, pi = pi_low, rho = rho)
}
ml_res <- data.frame()
for (iM in 1:4) {
  m <- multiLCA(dd, items, iT = K, id_high = if (iM > 1) "school" else NULL,
                iM = if (iM > 1) iM else NULL, incomplete = TRUE,
                verbose = FALSE, maxIter = 5000)
  mm <- lca_match(ref, to_fit(m))
  ml_res <- rbind(ml_res, data.frame(
    school_classes = iM,
    BIC_low = if (iM > 1) m$BIClow else m$BIC,
    BIC_high = if (iM > 1) m$BIChigh else NA,
    mad_vs_weighted = round(mm$mad, 4),
    shares = paste(sprintf("%.3f", mm$fit$pi), collapse = "/")))
  if (iM > 1) {
    cat(sprintf("\nC. %d school-level classes: school class sizes %s\n", iM,
                paste(sprintf("%.2f", m$vOmega), collapse = "/")))
    pg <- m$mPi[mm$perm, , drop = FALSE]
    rownames(pg) <- class_names
    print(round(pg, 3))
  }
}
cat("\nC. Multilevel LCA summary (1 = single level, package estimate):\n")
print(ml_res, row.names = FALSE)

write.csv(data.frame(class = class_names, weighted = ref$pi, unweighted = fu$pi),
          res_path("robust_unweighted_shares.csv"), row.names = FALSE)
write.csv(cbind(bw, unweighted = bu[, -1]), res_path("robust_unweighted_bch.csv"), row.names = FALSE)
write.csv(ml_res, res_path("robust_multilevel.csv"), row.names = FALSE)
