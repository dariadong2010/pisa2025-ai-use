# Correlates of class membership with missing correlate values multiply imputed (mice,
# 10 data sets) instead of complete cases. Steps 1-2 of the three-step method do not
# involve the correlates, so for each weight set (full + 80 BRR replicates) the LCA,
# assignment and D are computed once and step 3 is fitted on each imputed data set.
# Estimates are combined with Rubin's rules, using the BRR variance within each imputation.

source("code/common.R")
library(mice)

M_imp <- 10
ref <- load_reference()
d <- load_sample()
y <- as.matrix(d[, items])
oh_all <- one_hot(y, R)

# ---- Imputation ----
d$grade_below <- as.numeric(d$GRADE < 0)
d$eval_ai <- d$class_eval_ai - 1
d$private <- as.numeric(d$PRIVATESCH == 2)
d$community <- factor(ifelse(d$SC001Q01TA <= 3, "town", ifelse(d$SC001Q01TA == 4, "city", "large")),
                      levels = c("town", "city", "large"))
imp_dat <- data.frame(
  ESCS = d$ESCS, female = d$female, grade_below = d$grade_below, eval_ai = d$eval_ai,
  ai_school_avail = factor(d$ai_school_avail), ai_home_avail = factor(d$ai_home_avail),
  private = d$private, community = d$community,
  # auxiliary variables
  u1 = d$ST438Q01DA, u2 = d$ST438Q02DA, u3 = d$ST438Q03DA, u4 = d$ST438Q04DA,
  CURIO = d$CURIO, PERSEV = d$PERSEV, COGABIL = d$COGABIL, ENPROBS = d$ENPROBS,
  SCIE1 = d$PV1SCIE
)
cat("missing (%):\n"); print(round(100 * colMeans(is.na(imp_dat)), 1))
meth <- make.method(imp_dat)
meth[c("ESCS", "eval_ai", "u1", "u2", "u3", "u4", "CURIO", "PERSEV", "COGABIL", "ENPROBS")] <- "pmm"
imp <- mice(imp_dat, m = M_imp, method = meth, maxit = 10, seed = 20261101, printFlag = FALSE)
preds <- c("ESCS", "female", "grade_below", "eval_ai", "ai_school_avail", "ai_home_avail",
           "private", "city", "large_city")
X_list <- lapply(seq_len(M_imp), function(m) {
  cd <- complete(imp, m)
  X <- cbind(1, cd$ESCS, cd$female, cd$grade_below, cd$eval_ai,
             as.numeric(as.character(cd$ai_school_avail)), as.numeric(as.character(cd$ai_home_avail)),
             cd$private, as.numeric(cd$community == "city"), as.numeric(cd$community == "large"))
  colnames(X) <- c("(Intercept)", preds)
  X
})
P <- ncol(X_list[[1]])

# ---- Step 3 (as in 05_correlates_three_step.R) ----
mlogit_probs <- function(B, X) {
  eta <- X %*% t(B)
  eta <- eta - apply(eta, 1, max)
  p <- exp(eta)
  p / rowSums(p)
}
step3_fit <- function(X, Dw, w, start) {
  unpack <- function(b) rbind(0, matrix(b, K - 1, P))
  nll <- function(b) -sum(w * log(rowSums(mlogit_probs(unpack(b), X) * Dw)))
  gr <- function(b) {
    pr <- mlogit_probs(unpack(b), X)
    post <- pr * Dw / rowSums(pr * Dw)
    -as.vector((t(post - pr) %*% (X * w))[-1, ])
  }
  o <- optim(start, nll, gr, method = "BFGS", control = list(maxit = 2000, reltol = 1e-14))
  list(B = unpack(o$par), par = o$par, conv = o$convergence)
}

# One weight set -> list over imputations of fitted step-3 models
run_weightset <- function(wraw, start_b) {
  w <- wraw / mean(wraw)
  pat <- collapse_patterns(y, w)
  f <- lca_em(one_hot(pat$y, R), pat$w, K, R, start = ref, tol = 1e-11, maxit = 20000)
  f <- lca_match(ref, f)$fit
  post <- lca_posterior(f, oh_all)
  modal <- max.col(post, ties.method = "first")
  D <- t(sapply(seq_len(K), function(k)
    sapply(seq_len(K), function(s) sum(w * post[, k] * (modal == s)) / sum(w * post[, k]))))
  Dw <- t(D[, modal])
  lapply(seq_len(M_imp), function(m) step3_fit(X_list[[m]], Dw, w, start_b[[m]]))
}

est <- run_weightset(d$W_FSTUWT, rep(list(rep(0, (K - 1) * P)), M_imp))
starts <- lapply(est, `[[`, "par")
reps <- vector("list", 80)
for (r in 1:80) {
  reps[[r]] <- run_weightset(d[[paste0("W_FSTURWT", r)]], starts)
  if (r %% 20 == 0) cat("replicate", r, "\n")
}
cat("non-converged fits:", sum(unlist(lapply(c(list(est), reps), function(z) sapply(z, `[[`, "conv"))) != 0), "\n")

# ---- Rubin's rules with BRR variances, contrasts against a reference class ----
contrast_mi <- function(ref_class) {
  g <- function(B) as.vector(B - matrix(B[ref_class, ], K, P, byrow = TRUE))
  q <- sapply(est, function(z) g(z$B))                                   # parameters x imputations
  U <- sapply(seq_len(M_imp), function(m) {
    qr <- sapply(reps, function(z) g(z[[m]]$B))
    rowSums((qr - q[, m])^2) / (80 * 0.25)
  })
  qbar <- rowMeans(q)
  se <- sqrt(rowMeans(U) + (1 + 1 / M_imp) * apply(q, 1, var))
  out <- data.frame(reference = class_names[ref_class], class = rep(class_names, P),
                    term = rep(c("(Intercept)", preds), each = K), logit = qbar, se = se)
  out <- out[out$class != class_names[ref_class], ]
  out$p <- 2 * pnorm(-abs(out$logit / out$se))
  out$OR <- exp(out$logit); out$OR_lo <- exp(out$logit - 1.96 * out$se)
  out$OR_hi <- exp(out$logit + 1.96 * out$se)
  out
}
tab <- rbind(contrast_mi(1), contrast_mi(3))
write.csv(tab, res_path("correlates_or_mi.csv"), row.names = FALSE)

# Comparison with the complete-case results
cc <- read.csv(res_path("correlates_or.csv"))
cmp <- merge(tab[, c("reference", "class", "term", "OR", "p")],
             cc[, c("reference", "class", "term", "OR", "p")],
             by = c("reference", "class", "term"), suffixes = c("_mi", "_cc"))
cmp <- cmp[cmp$term != "(Intercept)", ]
cmp$sig_change <- (cmp$p_mi < .05) != (cmp$p_cc < .05)
write.csv(cmp, res_path("correlates_mi_vs_complete_cases.csv"), row.names = FALSE)
cat("\nOdds ratios, multiple imputation vs complete cases (rows where significance at .05 differs):\n")
show <- cmp[cmp$sig_change, ]
show[, c("OR_mi", "OR_cc")] <- round(show[, c("OR_mi", "OR_cc")], 2)
show[, c("p_mi", "p_cc")] <- signif(show[, c("p_mi", "p_cc")], 2)
print(show, row.names = FALSE)
