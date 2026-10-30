# RQ3 and H4-H5: class differences in learning dispositions and achievement with the
# BCH method (Bakk & Vermunt, 2016), using the six-class solution.
#   - Most likely class W and classification-error matrix D[k, s] = P(W = s | X = k);
#     each student enters class k with weight w_i * (D^-1)[W_i, k].
#   - Unadjusted class means, and class means adjusted for ESCS, gender and grade
#     (covariates centred, so class coefficients are means at the sample average).
#   - Achievement: each of the 10 plausible values separately, combined with Rubin's rules.
#   - SEs: the LCA, D and BCH steps are repeated for the 80 BRR-Fay replicate weights.
#   - Cohen's d against non-users uses the pooled within-class SD.
# The procedure is run on the full sample and again without possibly careless respondents.
# The disposition composite on the scale used in the paper (composite SD units) and the
# computational problem-solving results for all students are in 08_outcomes_composite.R.

source("code/common.R")

ref <- load_reference()
d0 <- load_sample()
d0$grade_below <- as.numeric(d0$GRADE < 0)
d0$careless <- (!is.na(d0$ST438_TT) & d0$ST438_TT < 8) | (!is.na(d0$q_effort) & d0$q_effort <= 5) |
  (d0$curio_straight %in% 1 & d0$persev_straight %in% 1)
covs <- c("ESCS", "female", "grade_below")
disp <- c("CURIO", "PERSEV", "COGABIL", "ENPROBS")
pv <- function(dom) paste0("PV", 1:M, dom)

# Standardisation constants from the full sample with final weights
wz <- d0$W_FSTUWT
zstd <- function(x) {
  ok <- !is.na(x)
  m <- weighted.mean(x[ok], wz[ok])
  (x - m) / sqrt(weighted.mean((x[ok] - m)^2, wz[ok]))
}
zd <- sapply(disp, function(v) zstd(d0[[v]]))
d0$DISP_z <- ifelse(rowSums(!is.na(zd)) >= 3, rowMeans(zd, na.rm = TRUE), NA)
for (m in 1:M) {
  d0[[paste0("SCIE_z", m)]] <- zstd(d0[[paste0("PV", m, "SCIE")]])
  d0[[paste0("GAP_SCIE", m)]] <- d0$DISP_z - d0[[paste0("SCIE_z", m)]]
  tk <- d0$ldw_taker == 1
  cz <- rep(NA_real_, nrow(d0))
  cz[tk] <- {
    x <- d0[[paste0("PV", m, "CMPS")]][tk]; ww <- wz[tk]; mu <- weighted.mean(x, ww)
    (x - mu) / sqrt(weighted.mean((x - mu)^2, ww))
  }
  d0[[paste0("GAP_CMPS", m)]] <- d0$DISP_z - cz
}

# Outcomes: variable names per plausible value (single variables repeated) and subset
outcomes <- list(
  CURIO = list(v = rep("CURIO", M)), PERSEV = list(v = rep("PERSEV", M)),
  COGABIL = list(v = rep("COGABIL", M)), ENPROBS = list(v = rep("ENPROBS", M)),
  DISP_z = list(v = rep("DISP_z", M)),
  SCIE = list(v = pv("SCIE")), READ = list(v = pv("READ")), MATH = list(v = pv("MATH")),
  CMPS_takers = list(v = pv("CMPS"), takers = TRUE), CMPS_all = list(v = pv("CMPS")),
  SCIE_z = list(v = paste0("SCIE_z", 1:M)),
  GAP_SCIE = list(v = paste0("GAP_SCIE", 1:M)),
  GAP_CMPS_takers = list(v = paste0("GAP_CMPS", 1:M), takers = TRUE)
)

# BCH class means for one outcome; the adjusted version solves the weighted normal
# equations directly (rows of D^-1 sum to 1, so the covariate block needs no H)
bch_means <- function(z, w, H, C = NULL) {
  ok <- !is.na(z)
  if (!is.null(C)) ok <- ok & complete.cases(C)
  wh <- w[ok] * H[ok, , drop = FALSE]
  if (is.null(C)) return(colSums(wh * z[ok]) / colSums(wh))
  Cc <- sweep(C[ok, , drop = FALSE], 2, colSums(w[ok] * C[ok, , drop = FALSE]) / sum(w[ok]))
  A <- rbind(cbind(diag(colSums(wh)), t(wh) %*% Cc),
             cbind(t(Cc) %*% wh, t(Cc) %*% (w[ok] * Cc)))
  b <- c(colSums(wh * z[ok]), colSums(w[ok] * z[ok] * Cc))
  solve(A, b)[1:K]
}

bch_pooled_sd <- function(z, w, H) {
  ok <- !is.na(z)
  wh <- w[ok] * H[ok, , drop = FALSE]
  mu <- colSums(wh * z[ok]) / colSums(wh)
  v <- colSums(wh * outer(z[ok], mu, "-")^2) / colSums(wh)
  sqrt(sum(colSums(wh) / sum(wh) * v))
}

# One weight set -> array [outcome, type (unadj/adj), class, plausible value]
run_bch <- function(d, wraw) {
  w <- wraw / mean(wraw)
  y <- as.matrix(d[, items])
  pat <- collapse_patterns(y, w)
  f <- lca_em(one_hot(pat$y, R), pat$w, K, R, start = ref, tol = 1e-11, maxit = 20000)
  f <- lca_match(ref, f)$fit
  post <- lca_posterior(f, one_hot(y, R))
  modal <- max.col(post, ties.method = "first")
  D <- t(sapply(seq_len(K), function(k)
    sapply(seq_len(K), function(s) sum(w * post[, k] * (modal == s)) / sum(w * post[, k]))))
  H <- solve(D)[modal, ]
  C <- as.matrix(d[, covs])
  out <- array(NA, c(length(outcomes), 2, K, M),
               dimnames = list(names(outcomes), c("unadj", "adj"), class_names, NULL))
  for (o in names(outcomes)) {
    sub <- if (isTRUE(outcomes[[o]]$takers)) d$ldw_taker == 1 else rep(TRUE, nrow(d))
    for (m in 1:M) {
      if (m > 1 && length(unique(outcomes[[o]]$v)) == 1) { out[o, , , m] <- out[o, , , 1]; next }
      z <- ifelse(sub, d[[outcomes[[o]]$v[m]]], NA)
      out[o, "unadj", , m] <- bch_means(z, w, H)
      out[o, "adj", , m] <- bch_means(z, w, H, C)
    }
  }
  list(means = out, pi = f$pi, H = H, w = w)
}

# Rubin's rules with BRR variances for any linear function g of the class means
combine <- function(est, reps, g) {
  q <- matrix(apply(est$means, 4, g), ncol = M)                   # quantities x PV
  qr <- array(sapply(reps, function(r) apply(r$means, 4, g), simplify = "array"),
              c(nrow(q), M, length(reps)))
  U <- rowMeans(apply(sweep(qr, 1:2, q)^2, 1:2, sum) / (length(reps) * 0.25))
  list(est = rowMeans(q), se = sqrt(U + (1 + 1 / M) * apply(q, 1, var)))
}

analyse <- function(d, label) {
  cat("\n==========", label, ": n =", nrow(d), "==========\n")
  est <- run_bch(d, d$W_FSTUWT)
  reps <- vector("list", 80)
  for (r in 1:80) {
    reps[[r]] <- run_bch(d, d[[paste0("W_FSTURWT", r)]])
    if (r %% 20 == 0) cat("replicate", r, "\n")
  }
  mean_tab <- do.call(rbind, lapply(names(outcomes), function(o) do.call(rbind, lapply(c("unadj", "adj"), function(t) {
    cm <- combine(est, reps, function(a) a[o, t, ])
    data.frame(sample = label, outcome = o, type = t, class = class_names, mean = cm$est, se = cm$se)
  }))))
  # All pairwise differences, Holm-adjusted within outcome and type
  pairs <- t(combn(K, 2))
  diff_tab <- do.call(rbind, lapply(names(outcomes), function(o) do.call(rbind, lapply(c("unadj", "adj"), function(t) {
    cm <- combine(est, reps, function(a) a[o, t, pairs[, 2]] - a[o, t, pairs[, 1]])
    p <- 2 * pnorm(-abs(cm$est / cm$se))
    data.frame(sample = label, outcome = o, type = t,
               class_a = class_names[pairs[, 2]], class_b = class_names[pairs[, 1]],
               diff = cm$est, se = cm$se, p = p, p_holm = p.adjust(p, "holm"))
  }))))
  d_tab <- do.call(rbind, lapply(names(outcomes), function(o) {
    sub <- if (isTRUE(outcomes[[o]]$takers)) d$ldw_taker == 1 else rep(TRUE, nrow(d))
    sds <- sapply(1:M, function(m) bch_pooled_sd(ifelse(sub, d[[outcomes[[o]]$v[m]]], NA), est$w, est$H))
    dd <- sapply(1:M, function(m) (est$means[o, "unadj", , m] - est$means[o, "unadj", 1, m]) / sds[m])
    data.frame(sample = label, outcome = o, class = class_names, cohen_d = rowMeans(dd))
  }))
  list(means = mean_tab, diffs = diff_tab, d = d_tab)
}

res_main <- analyse(d0, "main")
res_care <- analyse(d0[!d0$careless, ], "careless_excluded")
means <- rbind(res_main$means, res_care$means)
diffs <- rbind(res_main$diffs, res_care$diffs)
dtab <- rbind(res_main$d, res_care$d)
write.csv(means, res_path("outcomes_class_means.csv"), row.names = FALSE)
write.csv(diffs, res_path("outcomes_pairwise.csv"), row.names = FALSE)
write.csv(dtab, res_path("outcomes_cohen_d.csv"), row.names = FALSE)

for (s in c("main", "careless_excluded")) {
  cat("\n#####", s, "#####\n")
  for (t in c("unadj", "adj")) {
    cat("\nClass means (SE),", t, "\n")
    mt <- means[means$sample == s & means$type == t, ]
    tab <- sapply(names(outcomes), function(o) {
      x <- mt[mt$outcome == o, ]
      dig <- if (o %in% c("SCIE", "READ", "MATH", "CMPS_takers", "CMPS_all")) 0 else 2
      sprintf(paste0("%.", dig, "f (%.", dig, "f)"), x$mean, x$se)
    })
    rownames(tab) <- class_names
    print(noquote(t(tab)))
  }
  cat("\nCohen's d vs non-users (unadjusted):\n")
  dt <- dtab[dtab$sample == s, ]
  print(round(xtabs(cohen_d ~ outcome + class, dt)[names(outcomes)[1:10], class_names], 2))
}
