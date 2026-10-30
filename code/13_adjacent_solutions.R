# Key BCH class means under the five-, six- and seven-class solutions (final weight,
# point estimates): profiles, disposition composite (composite SD units), the composite
# adjusted for the ICT-item ERS index, ICT-item ERS and science (mean of 10 PVs).

source("code/common.R")

fits <- readRDS(res_path("lca_fits.rds"))$fits
d <- load_sample(response_style = TRUE)
y <- as.matrix(d[, items]); w0 <- d$W_FSTUWT
wmean <- function(x, w) { ok <- !is.na(x); sum(w[ok] * x[ok]) / sum(w[ok]) }
wsd <- function(x, w) { ok <- !is.na(x); m <- wmean(x, w); sqrt(sum(w[ok] * (x[ok] - m)^2) / sum(w[ok])) }
zw <- function(x) (x - wmean(x, w0)) / wsd(x, w0)
disp <- c("CURIO", "PERSEV", "COGABIL", "ENPROBS")
Z <- sapply(disp, function(v) zw(d[[v]]))
comp <- ifelse(rowSums(!is.na(Z)) >= 3, rowMeans(Z, na.rm = TRUE), NA)
sdc <- wsd(comp, w0); d$DISP <- comp / sdc
Rz <- sapply(seq_along(disp), function(j) residuals(lm(Z[, j] ~ d$ERS_ict + d$ARS_bal, weights = w0, na.action = na.exclude)))
d$DISP_adj_ict <- ifelse(rowSums(!is.na(Rz)) >= 3, rowMeans(Rz, na.rm = TRUE), NA) / sdc
d$SCIE <- rowMeans(sapply(1:M, function(m) d[[paste0("PV", m, "SCIE")]]))

for (Kc in 5:7) {
  f <- fits[[Kc]]
  post <- lca_posterior(f, one_hot(y, R)); modal <- max.col(post, ties.method = "first")
  w <- w0 / mean(w0)
  D <- t(sapply(seq_len(Kc), function(k) sapply(seq_len(Kc), function(s) sum(w * post[, k] * (modal == s)) / sum(w * post[, k]))))
  H <- solve(D)[modal, ]
  bm <- function(z) { ok <- !is.na(z); colSums(w[ok] * H[ok, ] * z[ok]) / colSums(w[ok] * H[ok, ]) }
  lev <- sapply(f$rho, function(r) r %*% seq_len(R))
  tab <- data.frame(class = seq_len(Kc), share = round(f$pi, 3), round(lev, 2),
                    DISP = round(bm(d$DISP), 2), DISP_adj_ict = round(bm(d$DISP_adj_ict), 2),
                    ERS_ict = round(bm(d$ERS_ict), 3), SCIE = round(bm(d$SCIE)))
  names(tab)[3:6] <- c("summ", "res", "draft", "learn")
  cat("\nK =", Kc, "\n"); print(tab, row.names = FALSE)
  write.csv(tab, res_path(sprintf("adjacent_K%d.csv", Kc)), row.names = FALSE)
}
