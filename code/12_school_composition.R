# Profile composition across schools (Figure S2).
#   1. Multilevel LCA with three school-level classes (multilevLCA, unweighted): sizes of the
#      school classes and the profile composition within each.
#   2. Estimated profile shares per school from the weighted single-level solution
#      (weighted mean of the students' posterior class probabilities within each school).

source("code/common.R")
library(multilevLCA)

set.seed(20261104)
ref <- load_reference()
d <- load_sample()
y <- as.matrix(d[, items])

# 1. Multilevel model
dd <- data.frame(y - 1, school = d$CNTSCHID)
names(dd)[1:4] <- items
m <- multiLCA(dd, items, iT = K, id_high = "school", iM = 3, incomplete = TRUE,
              verbose = FALSE, maxIter = 5000)
rho <- lapply(seq_along(items), function(j) t(m$mPhi[(j - 1) * R + seq_len(R), , drop = FALSE]))
mm <- lca_match(ref, list(K = K, pi = as.vector(m$mPi %*% m$vOmega), rho = rho))
comp <- m$mPi[mm$perm, , drop = FALSE]
ml <- data.frame(school_class = rep(seq_len(3), each = K),
                 school_class_share = rep(as.vector(m$vOmega), each = K),
                 class = rep(class_names, 3), share = as.vector(comp))
write.csv(ml, res_path("school_classes_composition.csv"), row.names = FALSE)
print(ml)

# 2. Per-school shares from the weighted solution
post <- lca_posterior(ref, one_hot(y, R))
w <- d$W_FSTUWT
sch <- do.call(rbind, lapply(split(seq_len(nrow(d)), d$CNTSCHID), function(i) {
  data.frame(CNTSCHID = d$CNTSCHID[i[1]], n = length(i),
             t(colSums(w[i] * post[i, , drop = FALSE]) / sum(w[i])))
}))
names(sch)[3:(2 + K)] <- class_names
write.csv(sch[, -1], res_path("school_profile_shares.csv"), row.names = FALSE)
cat("schools:", nrow(sch), " median students per school:", median(sch$n), "\n")
print(round(apply(sch[, class_names], 2, quantile, c(.1, .5, .9)), 3))
