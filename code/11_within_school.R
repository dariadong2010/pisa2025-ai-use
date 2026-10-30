# Class differences within schools: BCH-weighted regression of the outcome on class
# dummies with school fixed effects (within transformation). Rows are student x class with
# weight w_i * (D^-1)[W_i, k]; because each student's BCH weights sum to 1, school means of
# the outcome are plain weighted school means, and school means of a class dummy are the
# school's estimated class shares. Coefficients are within-school differences from non-users.
# SEs: LCA, D and BCH repeated for the 80 BRR-Fay replicates; plausible values by Rubin's rules.

source("code/common.R")

ref <- load_reference()
d <- load_sample()
y <- as.matrix(d[, items])
oh_all <- one_hot(y, R)
school <- as.integer(factor(d$CNTSCHID))

wz <- d$W_FSTUWT
zstd <- function(x, sub = rep(TRUE, length(x))) {
  ok <- !is.na(x) & sub; m <- weighted.mean(x[ok], wz[ok])
  out <- (x - m) / sqrt(weighted.mean((x[ok] - m)^2, wz[ok]))
  out[!sub] <- NA
  out
}
zd <- sapply(c("CURIO", "PERSEV", "COGABIL", "ENPROBS"), function(v) zstd(d[[v]]))
d$DISP_z <- ifelse(rowSums(!is.na(zd)) >= 3, rowMeans(zd, na.rm = TRUE), NA)
tk <- d$ldw_taker == 1
for (m in 1:M) {
  d[[paste0("GAP_SCIE", m)]] <- d$DISP_z - zstd(d[[paste0("PV", m, "SCIE")]])
  d[[paste0("CMPS_tk", m)]] <- ifelse(tk, d[[paste0("PV", m, "CMPS")]], NA)
}
outcomes <- list(DISP_z = rep("DISP_z", M), SCIE = paste0("PV", 1:M, "SCIE"),
                 GAP_SCIE = paste0("GAP_SCIE", 1:M), CMPS_takers = paste0("CMPS_tk", 1:M))

# Within-school BCH regression: K - 1 differences from class 1
within_fe <- function(z, w, H) {
  ok <- !is.na(z)
  z <- z[ok]; w <- w[ok]; H <- H[ok, , drop = FALSE]; s <- school[ok]
  sw <- tapply(w, s, sum)
  zbar <- as.numeric((tapply(w * z, s, sum) / sw)[as.character(s)])
  share <- apply(w * H[, -1, drop = FALSE], 2, function(v) (tapply(v, s, sum) / sw)[as.character(s)])
  zt <- z - zbar
  A <- matrix(0, K - 1, K - 1); b <- rep(0, K - 1)
  for (k in 1:K) {
    e <- matrix(0, length(z), K - 1)
    if (k > 1) e[, k - 1] <- 1
    Xt <- e - share
    om <- w * H[, k]
    A <- A + t(Xt) %*% (om * Xt)
    b <- b + colSums(om * zt * Xt)
  }
  solve(A, b)
}

run_ws <- function(wraw) {
  w <- wraw / mean(wraw)
  pat <- collapse_patterns(y, w)
  f <- lca_em(one_hot(pat$y, R), pat$w, K, R, start = ref, tol = 1e-11, maxit = 20000)
  f <- lca_match(ref, f)$fit
  post <- lca_posterior(f, oh_all)
  modal <- max.col(post, ties.method = "first")
  D <- t(sapply(seq_len(K), function(k)
    sapply(seq_len(K), function(s) sum(w * post[, k] * (modal == s)) / sum(w * post[, k]))))
  H <- solve(D)[modal, ]
  sapply(names(outcomes), function(o)
    sapply(1:M, function(m) {
      if (m > 1 && length(unique(outcomes[[o]])) == 1) return(rep(NA, K - 1))
      within_fe(d[[outcomes[[o]][m]]], w, H)
    }), simplify = "array")                       # (K - 1) x PV x outcome
}

est <- run_ws(d$W_FSTUWT)
reps <- vector("list", 80)
for (r in 1:80) {
  reps[[r]] <- run_ws(d[[paste0("W_FSTURWT", r)]])
  if (r %% 20 == 0) cat("replicate", r, "\n")
}

# Each class vs non-users, plus selected contrasts between user classes
L <- rbind(diag(K - 1),
           "moderate all-purpose vs learning-focused" = c(0, -1, 1, 0, 0),
           "daily all-purpose vs learning-focused" = c(0, -1, 0, 0, 1),
           "daily all-purpose vs moderate all-purpose" = c(0, 0, -1, 0, 1),
           "daily all-purpose vs frequent all-purpose" = c(0, 0, 0, -1, 1))
rownames(L)[1:(K - 1)] <- paste(class_names[-1], "vs non-user")
res <- do.call(rbind, lapply(names(outcomes), function(o) {
  Mo <- if (length(unique(outcomes[[o]])) == 1) 1 else M
  q <- matrix(L %*% est[, 1:Mo, o, drop = FALSE][, , 1], nrow(L))
  U <- rowMeans(matrix(sapply(1:Mo, function(m) {
    qr <- sapply(reps, function(z) L %*% z[, m, o])
    rowSums((qr - q[, m])^2) / (80 * 0.25)
  }), nrow(L)))
  B <- if (Mo > 1) apply(q, 1, var) else 0
  se <- sqrt(U + (1 + 1 / Mo) * B)
  data.frame(outcome = o, contrast = rownames(L), within_school_diff = rowMeans(q), se = se,
             p = 2 * pnorm(-abs(rowMeans(q) / se)))
}))
res$p_holm <- ave(res$p, res$outcome, FUN = function(p) p.adjust(p, "holm"))
write.csv(res, res_path("robust_within_school.csv"), row.names = FALSE)

show <- res
show$within_school_diff <- round(show$within_school_diff, 2); show$se <- round(show$se, 2)
show$p_holm <- signif(show$p_holm, 2); show$p <- signif(show$p, 2)
print(show, row.names = FALSE)
