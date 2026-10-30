# BCH class means on the scale used in the paper (80 BRR-Fay replicates, 10 plausible
# values with Rubin's rules), unadjusted and adjusted for ESCS, gender and grade.
#   - Disposition composite: mean of the four standardised indices divided by its weighted
#     SD, so 1 = one SD of the observed composite. Versions adjusted for response style are
#     the residuals of each standardised index on an ERS index and balanced acquiescence,
#     averaged and divided by the same SD.
#   - Difference between the composite and standardised science, and the composite at a
#     given science level (residual of the composite on science, per plausible value).
#   - Computational problem solving for all students and for students who took the domain.
#   - Share of the daily-vs-non-user composite difference that remains after adjustment.
# Run on the full sample and without students flagged by time/effort or by identical answers.

source("code/common.R")

ref <- load_reference()
d0 <- load_sample(response_style = TRUE)
d0$grade_below <- as.numeric(d0$GRADE < 0)
d0$flag_time_effort <- (!is.na(d0$ST438_TT) & d0$ST438_TT < 8) | (!is.na(d0$q_effort) & d0$q_effort <= 5)
d0$flag_pattern <- d0$curio_straight %in% 1 & d0$persev_straight %in% 1
covs <- c("ESCS", "female", "grade_below")
w0 <- d0$W_FSTUWT

wmean <- function(x, w) { ok <- !is.na(x); sum(w[ok] * x[ok]) / sum(w[ok]) }
wsd <- function(x, w) { ok <- !is.na(x); m <- wmean(x, w); sqrt(sum(w[ok] * (x[ok] - m)^2) / sum(w[ok])) }
zw <- function(x, w = w0) (x - wmean(x, w)) / wsd(x, w)

# ---- Outcomes (full sample, final weights) ----
disp <- c("CURIO", "PERSEV", "COGABIL", "ENPROBS")
Z <- sapply(disp, function(v) zw(d0[[v]]))
comp_raw <- ifelse(rowSums(!is.na(Z)) >= 3, rowMeans(Z, na.rm = TRUE), NA)
sd_comp <- wsd(comp_raw, w0)
write.csv(data.frame(composite_sd = sd_comp), res_path("composite_sd.csv"), row.names = FALSE)
d0$DISP <- comp_raw / sd_comp
for (v in c("ERS", "ERS_ict", "ERS_neu")) {
  Rz <- sapply(seq_along(disp), function(j)
    residuals(lm(Z[, j] ~ d0[[v]] + d0$ARS_bal, weights = w0, na.action = na.exclude)))
  d0[[paste0("DISP_adj_", v)]] <- ifelse(rowSums(!is.na(Rz)) >= 3, rowMeans(Rz, na.rm = TRUE), NA) / sd_comp
}
tk <- d0$ldw_taker == 1
for (m in 1:M) {
  sci <- zw(d0[[paste0("PV", m, "SCIE")]])
  d0[[paste0("SCIE_z", m)]] <- sci
  d0[[paste0("GAP", m)]] <- d0$DISP - sci
  d0[[paste0("GAP_adj_ict", m)]] <- d0$DISP_adj_ERS_ict - sci
  b <- coef(lm(d0$DISP ~ sci, weights = w0))[2]
  d0[[paste0("DISP_given_SCIE", m)]] <- d0$DISP - b * sci
  d0[[paste0("CMPS_tk", m)]] <- ifelse(tk, d0[[paste0("PV", m, "CMPS")]], NA)
}
pvs <- function(stem) paste0(stem, 1:M)
one <- function(v) rep(v, M)
outcomes <- list(
  ERS = one("ERS"), ERS_ict = one("ERS_ict"),
  DISP = one("DISP"), DISP_adj_ERS = one("DISP_adj_ERS"), DISP_adj_ERS_ict = one("DISP_adj_ERS_ict"),
  DISP_adj_ERS_neu = one("DISP_adj_ERS_neu"),
  SCIE = paste0("PV", 1:M, "SCIE"), SCIE_z = pvs("SCIE_z"),
  GAP = pvs("GAP"), GAP_adj_ict = pvs("GAP_adj_ict"), DISP_given_SCIE = pvs("DISP_given_SCIE"),
  CMPS_all = paste0("PV", 1:M, "CMPS"), CMPS_takers = pvs("CMPS_tk")
)

# ---- BCH (as in 07_outcomes_bch.R) ----
bch_means <- function(z, w, H, C = NULL) {
  ok <- !is.na(z)
  if (!is.null(C)) ok <- ok & complete.cases(C)
  wh <- w[ok] * H[ok, , drop = FALSE]
  if (is.null(C)) return(colSums(wh * z[ok]) / colSums(wh))
  Cc <- sweep(C[ok, , drop = FALSE], 2, colSums(w[ok] * C[ok, , drop = FALSE]) / sum(w[ok]))
  A <- rbind(cbind(diag(colSums(wh)), t(wh) %*% Cc), cbind(t(Cc) %*% wh, t(Cc) %*% (w[ok] * Cc)))
  solve(A, c(colSums(wh * z[ok]), colSums(w[ok] * z[ok] * Cc)))[1:K]
}
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
  for (o in names(outcomes)) for (m in 1:M) {
    if (m > 1 && length(unique(outcomes[[o]])) == 1) { out[o, , , m] <- out[o, , , 1]; next }
    z <- d[[outcomes[[o]][m]]]
    out[o, "unadj", , m] <- bch_means(z, w, H)
    out[o, "adj", , m] <- bch_means(z, w, H, C)
  }
  list(means = out, pi = f$pi)
}
combine <- function(est, reps, g) {
  q <- matrix(apply(est$means, 4, g), ncol = M)
  qr <- array(sapply(reps, function(r) apply(r$means, 4, g), simplify = "array"), c(nrow(q), M, length(reps)))
  U <- rowMeans(apply(sweep(qr, 1:2, q)^2, 1:2, sum) / (length(reps) * 0.25))
  list(est = rowMeans(q), se = sqrt(U + (1 + 1 / M) * apply(q, 1, var)))
}

analyse <- function(d, label) {
  cat("\n==========", label, ": n =", nrow(d), "==========\n")
  est <- run_bch(d, d$W_FSTUWT)
  reps <- lapply(1:80, function(r) {
    if (r %% 20 == 0) cat("replicate", r, "\n")
    run_bch(d, d[[paste0("W_FSTURWT", r)]])
  })
  means <- do.call(rbind, lapply(names(outcomes), function(o) do.call(rbind, lapply(c("unadj", "adj"), function(t) {
    cm <- combine(est, reps, function(a) a[o, t, ])
    data.frame(sample = label, outcome = o, type = t, class = class_names, mean = cm$est, se = cm$se)
  }))))
  pairs <- t(combn(K, 2))
  diffs <- do.call(rbind, lapply(names(outcomes), function(o) do.call(rbind, lapply(c("unadj", "adj"), function(t) {
    cm <- combine(est, reps, function(a) a[o, t, pairs[, 2]] - a[o, t, pairs[, 1]])
    p <- 2 * pnorm(-abs(cm$est / cm$se))
    data.frame(sample = label, outcome = o, type = t, class_a = class_names[pairs[, 2]],
               class_b = class_names[pairs[, 1]], diff = cm$est, se = cm$se, p = p, p_holm = p.adjust(p, "holm"))
  }))))
  rem <- do.call(rbind, lapply(c("DISP_adj_ERS", "DISP_adj_ERS_ict", "DISP_adj_ERS_neu"), function(a) {
    cm <- combine(est, reps, function(x) (x[a, "unadj", 6] - x[a, "unadj", 1]) / (x["DISP", "unadj", 6] - x["DISP", "unadj", 1]))
    data.frame(sample = label, adjusted = a, share_remaining = cm$est, se = cm$se)
  }))
  shares <- data.frame(sample = label, class = class_names, share = est$pi, n_students = nrow(d))
  list(means = means, diffs = diffs, rem = rem, shares = shares)
}

res <- list(analyse(d0, "main"),
            analyse(d0[!d0$flag_time_effort, ], "time_effort_excluded"),
            analyse(d0[!d0$flag_pattern, ], "pattern_excluded"))
for (k in c("means", "diffs", "rem", "shares"))
  write.csv(do.call(rbind, lapply(res, `[[`, k)), res_path(paste0("composite_", k, ".csv")), row.names = FALSE)

cat("\nweighted SD of the raw composite:", round(sd_comp, 3), "\n")
for (r in res) {
  s <- r$means$sample[1]
  cat("\n#####", s, "#####  class shares:", paste(sprintf("%.3f", r$shares$share), collapse = "/"), "\n")
  for (t in c("unadj", "adj")) {
    mt <- r$means[r$means$type == t, ]
    tab <- sapply(names(outcomes), function(o) {
      x <- mt[mt$outcome == o, ]; dg <- if (o %in% c("SCIE", "CMPS_all", "CMPS_takers")) 0 else 2
      sprintf(paste0("%.", dg, "f (%.", dg, "f)"), x$mean, x$se)
    })
    rownames(tab) <- class_names
    cat("\n", t, "\n"); print(noquote(t(tab)))
  }
  print(format(r$rem, digits = 3), row.names = FALSE)
  key <- r$diffs[r$diffs$type == "unadj" & r$diffs$class_a == "daily all-purpose" &
                 r$diffs$class_b %in% c("non-user", "frequent all-purpose", "learning-focused") &
                 r$diffs$outcome %in% c("DISP", "DISP_adj_ERS_ict", "SCIE", "CMPS_all", "GAP", "DISP_given_SCIE", "ERS_ict"), ]
  key$diff <- round(key$diff, 2); key$se <- round(key$se, 2); key$p_holm <- signif(key$p_holm, 2)
  print(key[, c("outcome", "class_b", "diff", "se", "p_holm")], row.names = FALSE)
}
