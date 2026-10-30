# Weighted latent class analysis for nominal indicators (pseudo-maximum likelihood).
#
# The likelihood sums w_i * log P(y_i) over students, with weights rescaled to sum to n
# (Patterson, Dayton & Graubard, 2002). Missing indicator values drop out of a student's
# likelihood (full-information ML under MAR). Students are collapsed into unique response
# patterns with summed weights, which makes the EM fast.

collapse_patterns <- function(y, w) {
  key <- apply(y, 1, function(r) paste(ifelse(is.na(r), "x", r), collapse = ""))
  agg <- tapply(w, key, sum)
  cnt <- tapply(rep(1, length(key)), key, sum)
  first <- match(names(agg), key)
  list(y = y[first, , drop = FALSE], w = as.numeric(agg), n = as.numeric(cnt),
       key = names(agg), row_key = key)
}

# One-hot arrays: list over items of m x R matrices (all zeros where the item is missing)
one_hot <- function(y, R) {
  lapply(seq_len(ncol(y)), function(j) {
    z <- matrix(0, nrow(y), R)
    ok <- !is.na(y[, j])
    z[cbind(which(ok), y[ok, j])] <- 1
    z
  })
}

# m x K matrix of log P(y | class k)
class_loglik <- function(oh, rho) {
  out <- 0
  for (j in seq_along(oh)) out <- out + oh[[j]] %*% t(log(rho[[j]]))
  out
}

lca_em <- function(oh, w, K, R, start = NULL, tol = 1e-9, maxit = 5000, floor = 1e-10) {
  J <- length(oh)
  if (is.null(start)) {
    pi_k <- as.numeric(rgamma(K, 1)); pi_k <- pi_k / sum(pi_k)
    rho <- lapply(seq_len(J), function(j) {
      r <- matrix(rgamma(K * R, 1), K, R); r / rowSums(r)
    })
  } else {
    pi_k <- start$pi; rho <- start$rho
  }
  ll_old <- -Inf
  for (it in seq_len(maxit)) {
    lp <- sweep(class_loglik(oh, rho), 2, log(pi_k), "+")
    mx <- apply(lp, 1, max)
    ll_i <- mx + log(rowSums(exp(lp - mx)))
    ll <- sum(w * ll_i)
    post <- exp(lp - ll_i)
    wp <- post * w
    pi_k <- colSums(wp) / sum(w)
    rho <- lapply(oh, function(z) {
      num <- t(wp) %*% z
      num <- pmax(num / rowSums(num), floor)
      num / rowSums(num)
    })
    if (ll - ll_old < tol * abs(ll)) break
    ll_old <- ll
  }
  list(K = K, pi = pi_k, rho = rho, ll = ll, iter = it, converged = it < maxit, post = post)
}

# Posterior class probabilities for any set of response patterns
lca_posterior <- function(fit, oh) {
  lp <- sweep(class_loglik(oh, fit$rho), 2, log(fit$pi), "+")
  mx <- apply(lp, 1, max)
  exp(lp - mx - log(rowSums(exp(lp - mx))))
}

# Many random starts; keeps the best solution and how often its log-likelihood recurred
lca_fit_multi <- function(oh, w, K, R, nstart = 100, ...) {
  fits <- lapply(seq_len(nstart), function(s) lca_em(oh, w, K, R, ...))
  ll <- vapply(fits, `[[`, numeric(1), "ll")
  best <- fits[[which.max(ll)]]
  best$ll_all <- ll
  best$best_replicated <- sum(abs(ll - max(ll)) < 1e-3)
  best
}

# Order classes by mean expected response across items (non-users first)
lca_order <- function(fit, R) {
  lev <- sapply(fit$rho, function(r) r %*% seq_len(R))
  o <- order(rowMeans(matrix(lev, nrow = fit$K)))
  fit$pi <- fit$pi[o]
  fit$rho <- lapply(fit$rho, function(r) r[o, , drop = FALSE])
  fit$post <- fit$post[, o, drop = FALSE]
  fit
}

lca_fit_stats <- function(fit, w, n, J, R) {
  K <- fit$K
  npar <- (K - 1) + K * J * (R - 1)
  ent <- if (K == 1) NA else {
    p <- pmax(fit$post, 1e-300)
    1 - sum(w * rowSums(-p * log(p))) / (n * log(K))
  }
  modal <- max.col(fit$post, ties.method = "first")
  avepp <- if (K == 1) NA else
    min(sapply(seq_len(K), function(k) weighted.mean(fit$post[modal == k, k], w[modal == k])))
  data.frame(K = K, LL = fit$ll, npar = npar,
             AIC = -2 * fit$ll + 2 * npar,
             BIC = -2 * fit$ll + npar * log(n),
             aBIC = -2 * fit$ll + npar * log((n + 2) / 24),
             CAIC = -2 * fit$ll + npar * (log(n) + 1),
             entropy = ent, min_class = min(fit$pi), min_avepp = avepp,
             best_ll_reps = fit$best_replicated, converged = fit$converged)
}

# Long table of item-response probabilities
lca_profile_table <- function(fit, items) {
  do.call(rbind, lapply(seq_along(items), function(j) {
    r <- fit$rho[[j]]
    data.frame(K = fit$K, class = rep(seq_len(fit$K), ncol(r)),
               class_share = rep(fit$pi, ncol(r)), item = items[j],
               category = rep(seq_len(ncol(r)), each = fit$K), prob = as.vector(r))
  }))
}

# Match the classes of fit b to fit a (same K) by brute-force permutation on the
# item-response probabilities; returns b reordered and the mean absolute difference
lca_match <- function(a, b) {
  K <- a$K
  perms <- function(v) if (length(v) == 1) list(v) else
    do.call(c, lapply(seq_along(v), function(i) lapply(perms(v[-i]), function(p) c(v[i], p))))
  va <- do.call(cbind, a$rho)
  vb <- do.call(cbind, b$rho)
  best <- NULL; best_d <- Inf
  for (p in perms(seq_len(K))) {
    d <- mean(abs(va - vb[p, , drop = FALSE]))
    if (d < best_d) { best_d <- d; best <- p }
  }
  b$pi <- b$pi[best]
  b$rho <- lapply(b$rho, function(r) r[best, , drop = FALSE])
  list(fit = b, mad = best_d, perm = best)
}
