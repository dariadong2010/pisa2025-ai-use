# (a) Continuous reference models (Table S8, lower panel): weighted regressions of the
#     disposition composite (composite SD units) and of science on use intensity and breadth,
#     with BRR SEs (and 10 PVs for science). Intensity = mean of the four purpose items
#     (1-5, not centred); intensity_sq = (intensity - 3)^2; breadth = number of purposes used
#     at least monthly. m1: intensity; m2: + intensity_sq + breadth; m3: m2 + ICT-item ERS
#     + balanced acquiescence. No background covariates.
# (b) Weighted regressions of the ERS indices (share of end-category answers) on age, grade,
#     gender and ESCS, with BRR SEs and 95% CIs.

source("code/common.R")

d <- load_sample(response_style = TRUE)
y <- as.matrix(d[, items]); w0 <- d$W_FSTUWT

wmean <- function(x, w) { ok <- !is.na(x); sum(w[ok] * x[ok]) / sum(w[ok]) }
wsd <- function(x, w) { ok <- !is.na(x); m <- wmean(x, w); sqrt(sum(w[ok] * (x[ok] - m)^2) / sum(w[ok])) }
disp <- c("CURIO", "PERSEV", "COGABIL", "ENPROBS")
Z <- sapply(disp, function(v) (d[[v]] - wmean(d[[v]], w0)) / wsd(d[[v]], w0))
comp <- ifelse(rowSums(!is.na(Z)) >= 3, rowMeans(Z, na.rm = TRUE), NA)
d$DISP <- comp / wsd(comp, w0)
d$intensity <- rowMeans(y, na.rm = TRUE)
d$breadth <- rowSums(y >= 3, na.rm = TRUE)
d$intensity_sq <- (d$intensity - 3)^2
d$grade_below <- as.numeric(d$GRADE < 0)

brr_combine <- function(fun, vals_fun, n_pv) {
  est_pv <- lapply(seq_len(n_pv), function(m) fun(w0, vals_fun(m)))
  est <- Reduce(`+`, est_pv) / n_pv
  U <- Reduce(`+`, lapply(seq_len(n_pv), function(m) {
    reps <- lapply(1:80, function(r) fun(d[[paste0("W_FSTURWT", r)]], vals_fun(m)))
    Reduce(`+`, lapply(reps, function(q) (q - est_pv[[m]])^2)) / 20
  })) / n_pv
  B <- if (n_pv > 1) Reduce(`+`, lapply(est_pv, function(q) (q - est)^2)) / (n_pv - 1) else 0
  list(est = est, se = sqrt(U + (1 + 1 / n_pv) * B))
}

# ---- (a) continuous models ----
models <- list(m1 = c("intensity"), m2 = c("intensity", "intensity_sq", "breadth"),
               m3 = c("intensity", "intensity_sq", "breadth", "ERS_ict", "ARS_bal"))
terms <- models$m3
fit_set <- function(wv, o) {
  sapply(models, function(xs) {
    dd <- data.frame(o = o, d[, xs, drop = FALSE], w = wv)
    cf <- coef(lm(as.formula(paste("o ~", paste(xs, collapse = " + "))), data = dd, weights = w))
    out <- setNames(rep(NA_real_, length(terms)), terms); out[xs] <- cf[xs]; out
  })
}
cont <- do.call(rbind, lapply(list(DISP = list(function(m) d$DISP, 1),
                                   SCIE = list(function(m) d[[paste0("PV", m, "SCIE")]], M)), function(z) {
  r <- brr_combine(fit_set, z[[1]], z[[2]])
  data.frame(term = rep(terms, length(models)), model = rep(names(models), each = length(terms)),
             est = as.vector(r$est), se = as.vector(r$se))
}))
cont$outcome <- rep(c("DISP", "SCIE"), each = length(terms) * length(models))
cont <- cont[!is.na(cont$est), c("outcome", "model", "term", "est", "se")]
cat("complete cases for the composite models:", sum(complete.cases(d[, c("DISP", terms)])), "\n")
write.csv(cont, res_path("continuous_models.csv"), row.names = FALSE)
print(format(cont, digits = 3), row.names = FALSE)

# ---- (b) correlates of the ERS indices ----
preds <- c("AGE", "grade_below", "female", "ESCS")
ers_fit <- function(wv, o) {
  dd <- data.frame(o = o, d[, preds], w = wv)
  coef(lm(as.formula(paste("o ~", paste(preds, collapse = " + "))), data = dd, weights = w))[preds]
}
ers_tab <- do.call(rbind, lapply(c("ERS_ict", "ERS_neu", "ERS"), function(v) {
  r <- brr_combine(ers_fit, function(m) d[[v]], 1)
  data.frame(index = v, term = preds, est = r$est, se = r$se,
             lo = r$est - 1.96 * r$se, hi = r$est + 1.96 * r$se, p = 2 * pnorm(-abs(r$est / r$se)))
}))
cat("\ncomplete cases for the ERS models:", sum(complete.cases(d[, c("ERS_ict", preds)])), "\n")
write.csv(ers_tab, res_path("ers_correlates.csv"), row.names = FALSE)
print(format(ers_tab, digits = 3), row.names = FALSE)
