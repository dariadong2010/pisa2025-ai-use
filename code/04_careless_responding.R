# Possibly careless respondents and the stability of the six-class solution.
# Flags: less than 6 or 8 seconds on the ST438 screen, self-rated questionnaire effort
# of 5 or lower (ST331Q03, 1-10), and identical answers to all items of both the ST301
# (curiosity) and ST307 (perseverance) scales. Each flagged group is dropped in turn
# and the model refitted from random starts.

source("code/common.R")

set.seed(20261030)
J <- length(items); n_start <- 100; tol <- 1e-11; maxit <- 20000
ref <- load_reference()

d <- load_sample()
y <- as.matrix(d[, items])
w <- d$W_FSTUWT / mean(d$W_FSTUWT)

flag <- data.frame(
  speed6 = !is.na(d$ST438_TT) & d$ST438_TT < 6,
  speed8 = !is.na(d$ST438_TT) & d$ST438_TT < 8,
  low_effort = !is.na(d$q_effort) & d$q_effort <= 5,
  straight = d$curio_straight %in% 1 & d$persev_straight %in% 1
)
flag$any <- flag$speed8 | flag$low_effort | flag$straight

# Numbers flagged and overlaps between the three criteria
counts <- data.frame(
  quantity = c("speed8", "low_effort", "straight", "any", "speed8 & low_effort",
               "speed8 & straight", "low_effort & straight", "all three"),
  n = c(sum(flag$speed8), sum(flag$low_effort), sum(flag$straight), sum(flag$any),
        sum(flag$speed8 & flag$low_effort), sum(flag$speed8 & flag$straight),
        sum(flag$low_effort & flag$straight), sum(flag$speed8 & flag$low_effort & flag$straight)))
print(counts, row.names = FALSE)
write.csv(counts, res_path("careless_flag_counts.csv"), row.names = FALSE)

# Share of each class (most likely assignment) that carries each flag
post <- lca_posterior(ref, one_hot(y, R))
modal <- max.col(post, ties.method = "first")
by_class <- sapply(flag, function(f) tapply(w * f, modal, sum) / tapply(w, modal, sum))
by_class <- rbind(by_class, all = sapply(flag, function(f) weighted.mean(f, w)))
cat("\nWeighted share flagged, by most likely class:\n")
print(round(by_class, 3))
write.csv(data.frame(class = rownames(by_class), round(by_class, 4)),
          res_path("careless_flags_by_class.csv"), row.names = FALSE)

# Refit after dropping each flagged group
res <- data.frame(); lev_out <- list()
for (s in names(flag)) {
  keep <- !flag[[s]]
  wk <- d$W_FSTUWT[keep] / mean(d$W_FSTUWT[keep])
  pk <- collapse_patterns(y[keep, ], wk)
  f <- lca_fit_multi(one_hot(pk$y, R), pk$w, K, R, nstart = n_start, tol = tol, maxit = maxit)
  m <- lca_match(ref, f)
  st <- lca_fit_stats(f, pk$w, sum(keep), J, R)
  res <- rbind(res, data.frame(excluded = s, n_dropped = sum(!keep), n = sum(keep),
                               mad_vs_full = m$mad, entropy = st$entropy,
                               best_ll_reps = f$best_replicated,
                               shares = paste(sprintf("%.3f", m$fit$pi), collapse = "/")))
  lev_out[[s]] <- round(sapply(m$fit$rho, function(r) r %*% seq_len(R)), 2)
}
cat("\nFull-sample shares:", paste(sprintf("%.3f", ref$pi), collapse = "/"), "\n")
print(res, row.names = FALSE)
cat("\nExpected levels per class (summarize, research, draft, learn) after dropping any flag:\n")
print(lev_out$any)
write.csv(res, res_path("careless_lca_refits.csv"), row.names = FALSE)
