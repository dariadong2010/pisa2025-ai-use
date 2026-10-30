# RQ2: correlates of class membership with the bias-adjusted three-step ML method
# (Vermunt, 2010), using the six-class solution.
#   step 1  weighted LCA (refitted from the full-sample solution for each weight set)
#   step 2  most likely class and the classification-error matrix D[k, s] = P(W = s | X = k)
#   step 3  weighted multinomial logit for P(X | z) with D as fixed misclassification
# Standard errors: the whole procedure is repeated for the 80 BRR-Fay replicate weights
# (Fay factor 0.5), so step-1 uncertainty is carried into the step-3 estimates.
# Three models: main (complete cases), possibly careless respondents excluded, and the
# main model with age and the ICT-item ERS index added.

source("code/common.R")

ref <- load_reference()
d_all <- load_sample(response_style = TRUE)
d_all$grade_below <- as.numeric(d_all$GRADE < 0)          # below the modal grade
d_all$eval_ai <- d_all$class_eval_ai - 1                  # 0 never .. 3 very often
d_all$private <- as.numeric(d_all$PRIVATESCH == 2)
d_all$city <- as.numeric(d_all$SC001Q01TA == 4)           # reference: town or smaller (1-3)
d_all$large_city <- as.numeric(d_all$SC001Q01TA >= 5)     # large city or megacity
d_all$age <- d_all$AGE - 15.5
wq <- d_all$W_FSTUWT
wz <- function(x) { ok <- !is.na(x); m <- sum(wq[ok] * x[ok]) / sum(wq[ok]); (x - m) / sqrt(sum(wq[ok] * (x[ok] - m)^2) / sum(wq[ok])) }
d_all$ERS_ict_z <- wz(d_all$ERS_ict)
d_all$careless <- (!is.na(d_all$ST438_TT) & d_all$ST438_TT < 8) | (!is.na(d_all$q_effort) & d_all$q_effort <= 5) |
  (d_all$curio_straight %in% 1 & d_all$persev_straight %in% 1)

base_preds <- c("ESCS", "female", "grade_below", "eval_ai", "ai_school_avail", "ai_home_avail",
                "private", "city", "large_city")
continuous <- c("ESCS", "eval_ai", "age", "ERS_ict_z")

mlogit_probs <- function(B, X) {
  eta <- X %*% t(B)
  eta <- eta - apply(eta, 1, max)
  p <- exp(eta)
  p / rowSums(p)
}

three_step_model <- function(d, preds, tag) {
  binary <- setdiff(preds, continuous)
  y <- as.matrix(d[, items])
  oh_all <- one_hot(y, R)
  cc <- complete.cases(d[, preds])
  X <- cbind(1, as.matrix(d[cc, preds]))
  colnames(X) <- c("(Intercept)", preds)
  P <- ncol(X)
  cat("\n==========", if (tag == "") "main" else tag, "==========\n")
  cat("LCA sample", nrow(d), "; complete correlates", sum(cc), sprintf("(%.1f%%)\n", 100 * mean(cc)))

  step3_fit <- function(Dw, w, start) {
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

  # Average marginal effects on class probabilities (percentage points):
  # derivative for continuous correlates, change from 0 to 1 for binary ones
  ame <- function(B, w) {
    sapply(preds, function(v) {
      j <- which(colnames(X) == v)
      if (v %in% binary) {
        X1 <- X; X1[, j] <- 1; X0 <- X; X0[, j] <- 0
        colSums(w * (mlogit_probs(B, X1) - mlogit_probs(B, X0))) / sum(w)
      } else {
        pr <- mlogit_probs(B, X)
        colSums(w * pr * sweep(matrix(B[, j], nrow(X), K, byrow = TRUE), 1, pr %*% B[, j])) / sum(w)
      }
    }) * 100
  }

  three_step <- function(wraw, start_b) {
    w <- wraw / mean(wraw)
    pat <- collapse_patterns(y, w)
    f <- lca_em(one_hot(pat$y, R), pat$w, K, R, start = ref, tol = 1e-11, maxit = 20000)
    f <- lca_match(ref, f)$fit
    post <- lca_posterior(f, oh_all)
    modal <- max.col(post, ties.method = "first")
    D <- t(sapply(seq_len(K), function(k)
      sapply(seq_len(K), function(s) sum(w * post[, k] * (modal == s)) / sum(w * post[, k]))))
    s3 <- step3_fit(t(D[, modal[cc]]), w[cc], start_b)
    list(B = s3$B, par = s3$par, conv = s3$conv, D = D, pi = f$pi, ame = ame(s3$B, w[cc]))
  }

  est <- three_step(d$W_FSTUWT, rep(0, (K - 1) * P))
  cat("Classification-error matrix D (rows: true class, columns: assigned):\n")
  print(round(est$D, 3))
  cat("Classification accuracy (weighted):", round(sum(diag(est$D) * est$pi), 3), "\n")
  reps <- vector("list", 80)
  for (r in 1:80) {
    reps[[r]] <- three_step(d[[paste0("W_FSTURWT", r)]], est$par)
    if (r %% 20 == 0) cat("replicate", r, "\n")
  }
  cat("non-converged replicates:", sum(sapply(reps, `[[`, "conv") != 0), "\n")
  brr_se <- function(full_val, rep_vals) sqrt(rowSums((rep_vals - full_val)^2) / (80 * 0.25))

  contrast_table <- function(ref_class) {
    b_full <- est$B - matrix(est$B[ref_class, ], K, P, byrow = TRUE)
    b_reps <- sapply(reps, function(z) as.vector(z$B - matrix(z$B[ref_class, ], K, P, byrow = TRUE)))
    out <- data.frame(reference = class_names[ref_class],
                      class = rep(class_names, P), term = rep(c("(Intercept)", preds), each = K),
                      logit = as.vector(b_full), se = brr_se(as.vector(b_full), b_reps))
    out <- out[out$class != class_names[ref_class], ]
    out$z <- out$logit / out$se
    out$p <- 2 * pnorm(-abs(out$z))
    out$OR <- exp(out$logit)
    out$OR_lo <- exp(out$logit - 1.96 * out$se)
    out$OR_hi <- exp(out$logit + 1.96 * out$se)
    out
  }
  tab_nonuser <- contrast_table(1)
  tab_learn <- contrast_table(3)
  write.csv(rbind(tab_nonuser, tab_learn), res_path(paste0("correlates_or", tag, ".csv")), row.names = FALSE)

  ame_reps <- sapply(reps, function(z) as.vector(z$ame))
  ame_tab <- data.frame(class = rep(class_names, length(preds)), term = rep(preds, each = K),
                        ame_pp = as.vector(est$ame), se = brr_se(as.vector(est$ame), ame_reps))
  ame_tab$p <- 2 * pnorm(-abs(ame_tab$ame_pp / ame_tab$se))
  write.csv(ame_tab, res_path(paste0("correlates_ame", tag, ".csv")), row.names = FALSE)

  show <- function(tab) {
    tab <- tab[tab$term != "(Intercept)", ]
    fmt <- sprintf("%.2f [%.2f, %.2f]%s", tab$OR, tab$OR_lo, tab$OR_hi,
                   ifelse(tab$p < .001, "***", ifelse(tab$p < .01, "**", ifelse(tab$p < .05, "*", ""))))
    print(noquote(matrix(fmt, nrow = length(preds), byrow = TRUE, dimnames = list(preds, unique(tab$class)))))
  }
  cat("\nOdds ratios vs non-users [95% CI]:\n"); show(tab_nonuser)
  cat("\nOdds ratios vs learning-focused [95% CI]:\n"); show(tab_learn)
  cat("\nAverage marginal effects (percentage points; SE):\n")
  print(noquote(t(matrix(sprintf("%.1f (%.1f)", ame_tab$ame_pp, ame_tab$se), nrow = K,
                         dimnames = list(class_names, preds)))))
}

three_step_model(d_all, base_preds, "")
three_step_model(d_all[!d_all$careless, ], base_preds, "_careless")
three_step_model(d_all, c(base_preds[1:3], "age", base_preds[4:9], "ERS_ict_z"), "_ers_age")
cat("\nMean age by grade below the modal grade (0/1):",
    round(tapply(d_all$AGE, d_all$grade_below, mean, na.rm = TRUE), 3), "\n")
