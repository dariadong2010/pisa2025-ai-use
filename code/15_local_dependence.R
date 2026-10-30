# Sensitivity of the six-class solution to the local dependence left in the item pairs with
# the largest bivariate residuals: draft written work x help me learn, and summarize a text x
# research a new topic.
#
# Model: within each class the two items of a relaxed pair follow a log-linear model with
# class-specific main effects and one class-invariant direct effect (Hagenaars, 1988); all
# other items stay locally independent. Two forms of the direct effect:
#   same    one parameter for giving the same answer to both items (diagonal cells)
#   linear  linear-by-linear association of the 1-5 scores
# Variants: a/b = linear for the draft pair / both pairs; c/d = same for the draft pair /
# both pairs. Estimation is weighted pseudo-ML by EM; the M-step for a pair alternates
# iterative proportional fitting of the class-specific margins with a Newton step for the
# direct effect.
# Reported: fit and bivariate residuals; class shares and profiles; agreement of the most
# likely classes with the main model; key BCH results with BRR SEs.

source("code/common.R")

ref <- load_reference()
d <- load_sample(response_style = TRUE)
y <- as.matrix(d[, items])
w0 <- d$W_FSTUWT
n <- nrow(d)
cat("students:", n, "; missing per item:", colSums(is.na(y)), "\n")
fits_li <- readRDS(res_path("lca_fits.rds"))$fits

# Outcomes on the composite scale used in the paper
wmean <- function(x, w) { ok <- !is.na(x); sum(w[ok] * x[ok]) / sum(w[ok]) }
wsd <- function(x, w) { ok <- !is.na(x); m <- wmean(x, w); sqrt(sum(w[ok] * (x[ok] - m)^2) / sum(w[ok])) }
disp <- c("CURIO", "PERSEV", "COGABIL", "ENPROBS")
Z <- sapply(disp, function(v) (d[[v]] - wmean(d[[v]], w0)) / wsd(d[[v]], w0))
comp_raw <- ifelse(rowSums(!is.na(Z)) >= 3, rowMeans(Z, na.rm = TRUE), NA)
d$DISP <- comp_raw / wsd(comp_raw, w0)
variants_ers <- c("ERS", "ERS_ict", "ERS_neu")
for (v in variants_ers) {
  Rz <- sapply(seq_along(disp), function(j)
    residuals(lm(Z[, j] ~ d[[v]] + d$ARS_bal, weights = w0, na.action = na.exclude)))
  d[[paste0("DISP_adj_", v)]] <- ifelse(rowSums(!is.na(Rz)) >= 3, rowMeans(Rz, na.rm = TRUE), NA) /
    wsd(comp_raw, w0)
}
d$SCIE <- rowMeans(sapply(1:M, function(m) d[[paste0("PV", m, "SCIE")]]))
outcomes <- c("DISP", paste0("DISP_adj_", variants_ers), variants_ers, "SCIE")

# ---- Latent class model with direct effects ----
setup <- function(y, pairs) {
  S <- list(m = nrow(y), oh = one_hot(y, R), singles = setdiff(seq_len(ncol(y)), unlist(pairs)), pair = list())
  for (p in seq_along(pairs)) {
    a <- y[, pairs[[p]][1]]; b <- y[, pairs[[p]][2]]
    both <- !is.na(a) & !is.na(b); aonly <- !is.na(a) & is.na(b); bonly <- is.na(a) & !is.na(b)
    cellmat <- matrix(0, nrow(y), R * R)
    cellmat[cbind(which(both), (a + (b - 1) * R)[both])] <- 1
    amat <- matrix(0, nrow(y), R); amat[cbind(which(aonly), a[aonly])] <- 1
    bmat <- matrix(0, nrow(y), R); bmat[cbind(which(bonly), b[bonly])] <- 1
    S$pair[[p]] <- list(cellmat = cellmat, amat = amat, bmat = bmat, none = is.na(a) & is.na(b))
  }
  S
}

pair_Q <- function(p) {                                   # K x R x R within-class joint
  Kp <- nrow(p$A); Q <- array(0, c(Kp, R, R))
  for (k in seq_len(Kp)) {
    e <- exp(outer(p$A[k, ], p$B[k, ], "+") + p$delta * p$G)
    Q[k, , ] <- e / sum(e)
  }
  Q
}

de_loglik <- function(par, S) {                           # m x K of log P(y | class)
  out <- 0
  for (j in S$singles) out <- out + S$oh[[j]] %*% t(log(par$rho[[j]]))
  for (p in seq_along(par$pairs)) {
    Q <- pair_Q(par$pairs[[p]]); s <- S$pair[[p]]
    Kp <- dim(Q)[1]
    out <- out + s$cellmat %*% t(log(matrix(Q, Kp))) +
      s$amat %*% t(log(apply(Q, c(1, 2), sum))) + s$bmat %*% t(log(apply(Q, c(1, 3), sum)))
  }
  out
}

pair_counts <- function(Q, s, wp) {                       # expected K x R x R cell counts
  Kp <- dim(Q)[1]
  N <- array(t(wp) %*% s$cellmat, c(Kp, R, R))
  mA <- apply(Q, c(1, 2), sum); mB <- apply(Q, c(1, 3), sum)
  ca <- t(wp) %*% s$amat; cb <- t(wp) %*% s$bmat; cn <- colSums(wp[s$none, , drop = FALSE])
  for (k in seq_len(Kp)) {
    N[k, , ] <- N[k, , ] + Q[k, , ] * (ca[k, ] / mA[k, ]) +
      sweep(Q[k, , ], 2, cb[k, ] / mB[k, ], "*") + cn[k] * Q[k, , ]
  }
  N
}

pair_mstep <- function(p, N, inner = 5) {
  Kp <- dim(N)[1]; tot <- apply(N, 1, sum)
  rN <- pmax(apply(N, c(1, 2), sum), 1e-10); cN <- pmax(apply(N, c(1, 3), sum), 1e-10)
  Gk <- rep(p$G, each = Kp)
  sN <- sum(N * Gk)
  for (it in seq_len(inner)) {
    for (s in 1:3) {
      F <- pair_Q(p) * tot
      p$A <- p$A + log(rN / apply(F, c(1, 2), sum))
      F <- pair_Q(p) * tot
      p$B <- p$B + log(cN / apply(F, c(1, 3), sum))
    }
    Q <- pair_Q(p)
    g <- sN - sum(Q * tot * Gk)
    v <- sum(tot * (apply(Q * Gk^2, 1, sum) - apply(Q * Gk, 1, sum)^2))
    p$delta <- p$delta + g / v
  }
  p$A <- p$A - rowMeans(p$A); p$B <- p$B - rowMeans(p$B)
  p
}

de_em <- function(S, w, par, tol = 1e-10, maxit = 20000) {
  ll_old <- -Inf
  for (it in seq_len(maxit)) {
    lp <- sweep(de_loglik(par, S), 2, log(par$pi), "+")
    mx <- apply(lp, 1, max)
    ll_i <- mx + log(rowSums(exp(lp - mx)))
    ll <- sum(w * ll_i)
    wp <- exp(lp - ll_i) * w
    par$pi <- colSums(wp) / sum(w)
    for (j in S$singles) {
      num <- t(wp) %*% S$oh[[j]]
      num <- pmax(num / rowSums(num), 1e-10)
      par$rho[[j]] <- num / rowSums(num)
    }
    for (p in seq_along(par$pairs))
      par$pairs[[p]] <- pair_mstep(par$pairs[[p]], pair_counts(pair_Q(par$pairs[[p]]), S$pair[[p]], wp))
    if (it > 2 && ll - ll_old < tol * abs(ll)) break
    ll_old <- ll
  }
  par$ll <- ll; par$iter <- it; par$converged <- it < maxit
  par
}

start_from <- function(f, pairs, G) {
  list(pi = f$pi, rho = f$rho, pairs = lapply(pairs, function(pr)
    list(items = pr, A = log(f$rho[[pr[1]]]), B = log(f$rho[[pr[2]]]), delta = 0, G = G)))
}
start_random <- function(Kc, pairs, G) {
  rr <- function() { r <- matrix(rgamma(Kc * R, 1), Kc, R); r / rowSums(r) }
  f <- list(pi = { p <- rgamma(Kc, 1); p / sum(p) }, rho = lapply(seq_along(items), function(j) rr()))
  start_from(f, pairs, G)
}

# Item margins per class (for profiles, ordering and matching)
marginals <- function(par) {
  rho <- par$rho
  for (p in par$pairs) {
    Q <- pair_Q(p)
    rho[[p$items[1]]] <- apply(Q, c(1, 2), sum); rho[[p$items[2]]] <- apply(Q, c(1, 3), sum)
  }
  rho
}
permute <- function(par, o) {
  par$pi <- par$pi[o]
  par$rho <- lapply(par$rho, function(r) r[o, , drop = FALSE])
  par$pairs <- lapply(par$pairs, function(p) { p$A <- p$A[o, , drop = FALSE]; p$B <- p$B[o, , drop = FALSE]; p })
  par
}
match_to <- function(target, par) {
  m <- lca_match(target, list(K = length(par$pi), pi = par$pi, rho = marginals(par)))
  out <- permute(par, m$perm); out$mad <- m$mad
  out
}
npar_de <- function(Kc, pairs) (Kc - 1) + (length(items) - 2 * length(pairs)) * Kc * (R - 1) +
  length(pairs) * (2 * Kc * (R - 1) + 1)
posterior_de <- function(par, S) {
  lp <- sweep(de_loglik(par, S), 2, log(par$pi), "+")
  mx <- apply(lp, 1, max)
  exp(lp - mx - log(rowSums(exp(lp - mx))))
}

# Bivariate residuals for all item pairs under a model given by its class-conditional joints
bvr_all <- function(pi_k, joint_fun, w) {
  res <- c()
  for (j in 1:3) for (l in (j + 1):4) {
    ok <- !is.na(y[, j]) & !is.na(y[, l])
    obs <- xtabs(w[ok] ~ factor(y[ok, j], 1:R) + factor(y[ok, l], 1:R))
    ex <- sum(obs) * Reduce(`+`, lapply(seq_along(pi_k), function(k) pi_k[k] * joint_fun(j, l, k)))
    res[paste(j, l, sep = "x")] <- sum((obs - ex)^2 / ex) / ((R - 1)^2)
  }
  res
}
joint_de <- function(par) {
  rho <- marginals(par); Qs <- lapply(par$pairs, pair_Q)
  function(j, l, k) {
    for (p in seq_along(par$pairs)) if (all(par$pairs[[p]]$items == c(j, l))) return(Qs[[p]][k, , ])
    outer(rho[[j]][k, ], rho[[l]][k, ])
  }
}

bch <- function(post, w) {
  Kc <- ncol(post); modal <- max.col(post, ties.method = "first")
  D <- t(sapply(seq_len(Kc), function(k)
    sapply(seq_len(Kc), function(s) sum(w * post[, k] * (modal == s)) / sum(w * post[, k]))))
  H <- solve(D)[modal, ]
  sapply(outcomes, function(o) {
    z <- d[[o]]; ok <- !is.na(z)
    colSums(w[ok] * H[ok, ] * z[ok]) / colSums(w[ok] * H[ok, ])
  })
}
key_stats <- function(e) {
  gap <- unname(e[6, "DISP"] - e[1, "DISP"])
  c(gap_disp = gap,
    setNames(sapply(variants_ers, function(v) unname(e[6, paste0("DISP_adj_", v)] - e[1, paste0("DISP_adj_", v)]) / gap),
             paste0("remain_", variants_ers)),
    sci_learn_vs_non = unname(e[3, "SCIE"] - e[1, "SCIE"]), sci_daily_vs_non = unname(e[6, "SCIE"] - e[1, "SCIE"]))
}

# ---- Main model (local independence) ----
w <- w0 / mean(w0)
pat_all <- collapse_patterns(y, w)
post_ref <- lca_posterior(ref, one_hot(y, R))
modal_ref <- max.col(post_ref, ties.method = "first")
bvr_ref <- bvr_all(ref$pi, function(j, l, k) outer(ref$rho[[j]][k, ], ref$rho[[l]][k, ]), w)
est_ref <- bch(post_ref, w)
cat("\nMain model: LL", round(ref$ll, 2), "; bivariate residuals", round(bvr_ref, 1), "\n")
cat("unweighted counts by most likely class:", table(factor(modal_ref, 1:K)), "\n")

pairs_draft <- list(c(3, 4)); pairs_both <- list(c(3, 4), c(1, 2))
specs <- list(
  linear = list(G = outer(1:R, 1:R), variants = list(a = pairs_draft, b = pairs_both)),
  same = list(G = diag(R), variants = list(c = pairs_draft, d = pairs_both)))
fits <- list(); summary_rows <- list(); profile_rows <- list(); agree_rows <- list(); key_rows <- list()
bvr_rows <- list(data.frame(variant = "main", t(bvr_ref)))

for (sp in names(specs)) {
  set.seed(20261029)
  G <- specs[[sp]]$G
  for (vn in names(specs[[sp]]$variants)) {
    pairs <- specs[[sp]]$variants[[vn]]
    S_pat <- setup(pat_all$y, pairs)
    for (Kc in 5:7) {
      starts <- c(list(start_from(fits_li[[Kc]], pairs, G)),
                  lapply(1:(if (Kc == 6) 30 else 15), function(s) start_random(Kc, pairs, G)))
      cand <- lapply(starts, function(st) de_em(S_pat, pat_all$w, st))
      lls <- sapply(cand, `[[`, "ll")
      best <- cand[[which.max(lls)]]
      best$ll_reps <- sum(abs(lls - max(lls)) < 1e-3)
      lev <- sapply(marginals(best), function(r) r %*% 1:R)
      best <- permute(best, order(rowMeans(matrix(lev, nrow = Kc))))
      if (Kc == 6) best <- match_to(ref, best)
      np <- npar_de(Kc, pairs)
      li <- fits_li[[Kc]]
      np_li <- (Kc - 1) + Kc * length(items) * (R - 1)
      summary_rows[[length(summary_rows) + 1]] <- data.frame(
        variant = vn, direct_effect = sp, K = Kc, LL = best$ll, npar = np, BIC = -2 * best$ll + np * log(n),
        LL_li = li$ll, BIC_li = -2 * li$ll + np_li * log(n),
        delta = paste(round(sapply(best$pairs, `[[`, "delta"), 3), collapse = "; "),
        best_ll_reps = best$ll_reps, li_start_reached_best = abs(lls[1] - best$ll) < 1e-3,
        converged = best$converged)
      fits[[paste0(vn, "_K", Kc)]] <- best
      cat(sprintf("\nvariant %s (%s), K = %d: LL %.2f (main %.2f), BIC %.1f (main %.1f), direct effect %s, best LL reached %d times\n",
                  vn, sp, Kc, best$ll, li$ll, -2 * best$ll + np * log(n), -2 * li$ll + np_li * log(n),
                  paste(round(sapply(best$pairs, `[[`, "delta"), 3), collapse = "; "), best$ll_reps))
    }

    f6 <- fits[[paste0(vn, "_K6")]]
    post <- posterior_de(f6, setup(y, pairs))
    modal <- max.col(post, ties.method = "first")
    bv <- bvr_all(f6$pi, joint_de(f6), w)
    bvr_rows[[vn]] <- data.frame(variant = vn, t(bv))
    cat("bivariate residuals:", round(bv, 1), "; mean abs. difference in item probabilities vs main:", round(f6$mad, 4), "\n")
    lev <- sapply(marginals(f6), function(r) r %*% 1:R)
    lev_ref <- sapply(ref$rho, function(r) r %*% 1:R)
    profile_rows[[vn]] <- data.frame(variant = vn, class = class_names, share = f6$pi, share_ref = ref$pi,
                                     matrix(lev, K, dimnames = list(NULL, c("summ", "res", "draft", "learn"))),
                                     matrix(lev_ref, K, dimnames = list(NULL, paste0(c("summ", "res", "draft", "learn"), "_ref"))),
                                     n_modal = as.vector(table(factor(modal, 1:K))),
                                     n_modal_ref = as.vector(table(factor(modal_ref, 1:K))))
    same <- modal == modal_ref
    agree_rows[[vn]] <- data.frame(variant = vn, class = c(class_names, "all"),
                                   same_class_weighted = c(sapply(1:K, function(k) sum(w[modal_ref == k & same]) / sum(w[modal_ref == k])),
                                                           sum(w[same]) / sum(w)))
    cat("Weighted cross-tabulation of most likely classes (rows main, columns relaxed), row %:\n")
    print(round(100 * prop.table(xtabs(w ~ factor(modal_ref, 1:K, labels = class_names) + factor(modal, 1:K)), 1), 1))
    fits[[paste0(vn, "_bch")]] <- bch(post, w)
  }
}

# ---- BRR SEs for the key statistics of each six-class variant ----
key_main <- key_stats(est_ref)
se_rows <- list(data.frame(variant = "main", stat = names(key_main), est = key_main, se = NA))
for (sp in names(specs)) for (vn in names(specs[[sp]]$variants)) {
  pairs <- specs[[sp]]$variants[[vn]]
  f6 <- fits[[paste0(vn, "_K6")]]
  S_all <- setup(y, pairs)
  est_k <- key_stats(fits[[paste0(vn, "_bch")]])
  reps <- sapply(1:80, function(r) {
    wr <- d[[paste0("W_FSTURWT", r)]]; wr <- wr / mean(wr)
    pr <- collapse_patterns(y, wr)
    fr <- de_em(setup(pr$y, pairs), pr$w, f6, tol = 1e-11)
    fr <- match_to(list(K = K, pi = f6$pi, rho = marginals(f6)), fr)
    if (r %% 20 == 0) cat("variant", vn, "replicate", r, "\n")
    key_stats(bch(posterior_de(fr, S_all), wr))
  })
  se_rows[[vn]] <- data.frame(variant = vn, stat = names(est_k), est = est_k,
                              se = sqrt(rowSums((reps - est_k)^2) / (80 * 0.25)))
}

out <- list(fit = do.call(rbind, summary_rows), bvr = do.call(rbind, bvr_rows),
            profiles = do.call(rbind, profile_rows), agreement = do.call(rbind, agree_rows),
            key = do.call(rbind, se_rows))
for (k in names(out)) write.csv(out[[k]], res_path(paste0("local_dependence_", k, ".csv")), row.names = FALSE)
for (k in names(out)) { cat("\n", k, ":\n"); print(format(out[[k]], digits = 4), row.names = FALSE) }
