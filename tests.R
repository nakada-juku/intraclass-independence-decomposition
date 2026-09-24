# =====================================================================
#  Checks. Run with:  Rscript tests.R
#
#  Three kinds of check:
#    (a) the values printed in the paper are reproduced;
#    (b) the fits satisfy the relations they must satisfy, on random tables;
#    (c) the public entry point accepts every documented form of input.
# =====================================================================

source("R/intraclass.R")
source("R/populations.R")

fails <- 0L
ok <- function(label, passed, detail = "") {
  cat(sprintf("%-58s %s%s\n", label, if (passed) "ok" else "FAIL",
              if (nzchar(detail)) paste0("  ", detail) else ""))
  if (!passed) fails <<- fails + 1L
}

## (a) the paper's numbers ---------------------------------------------
DATA <- list(
  "Ishii (1960)" = list(y = c(4, 17, 8, 21, 20, 6), r = 3,
                        G2 = c(ZC = 0.190, LL = 0.049, I = 0.243), theta = 1.112),
  "Data 1"       = list(y = c(14, 85, 194, 65, 81, 272, 136, 263, 318, 42), r = 4,
                        G2 = c(ZC = 0.670, LL = 31.275, I = 31.836), theta = NA),
  "Data 2"       = list(y = c(45, 85, 74, 33, 65, 93, 69, 62, 52, 24), r = 4,
                        G2 = c(ZC = 9.121, LL = 9.717, I = 18.904), theta = NA))
cat("-- (a) the values printed in the paper --\n")
for (nm in names(DATA)) {
  d <- DATA[[nm]]; f <- decompose(d$y, d$r)
  for (m in c("ZC", "LL", "I"))
    ok(sprintf("%s: G2(H_%s) = %.3f", nm, m, d$G2[[m]]),
       abs(f$G2[[m]] - d$G2[[m]]) < 0.005, sprintf("got %.4f", f$G2[[m]]))
  if (!is.na(d$theta))
    ok(sprintf("%s: theta = %.3f", nm, d$theta), abs(f$theta - d$theta) < 0.001)
}
ok("Ishii: fitted frequencies match the table in the paper",
   max(abs(decompose(DATA[[1]]$y, 3)$fits$ZC$m -
           c(3.594, 16.916, 8.919, 20.995, 20.091, 5.484))) < 0.001)

## (b) relations that must hold ----------------------------------------
cat("\n-- (b) relations that must hold, on 270 random tables --\n")
set.seed(20260922)
bad_nested <- bad_cov <- bad_sum <- bad_fit <- 0L; N <- 0L
maxgrad <- 0; worst_sum <- 0; worst_neg <- 0
for (r in c(3, 4, 5)) for (nn in c(40, 100, 400))
  for (sh in names(A_VECTORS)) for (cs in list(c(0, 0), c(0.3, 0.8))) {
    u <- seq_len(r); p <- pop_pi(A_VECTORS[[sh]](r), cs[1], cs[2], u)
    for (b in 1:5) {
      N <- N + 1L
      y <- as.vector(rmultinom(1, nn, p))
      f <- decompose(y, r, u)
      if (is.null(f)) { bad_fit <- bad_fit + 1L; next }
      # H_I is nested in H_ZC and in H_LL
      tol <- 1e-6 * max(1, f$G2[["I"]])
      if (f$G2[["ZC"]] > f$G2[["I"]] + tol || f$G2[["LL"]] > f$G2[["I"]] + tol)
        bad_nested <- bad_nested + 1L
      # the ZC fit satisfies its constraint
      if (abs(cov_UV(f$fits$ZC$pi, r, u)) > 1e-5) bad_cov <- bad_cov + 1L
      # every fit is a probability vector. The constrained ZC fit meets its
      # three linear constraints to the documented tolerance of 1e-6; the
      # attained accuracy is reported below.
      for (fit in f$fits) {
        worst_sum <<- max(worst_sum, abs(sum(fit$pi) - 1))
        worst_neg <<- min(worst_neg, min(fit$pi))
        if (abs(sum(fit$pi) - 1) > 1e-6 || any(fit$pi < -1e-12)) bad_sum <- bad_sum + 1L
      }
      maxgrad <- max(maxgrad, f$fits$LL$gradient_max)
    }
  }
ok(sprintf("no failed fit out of %d", N), bad_fit == 0L, sprintf("%d failed", bad_fit))
ok("G2(H_I) >= G2(H_ZC), G2(H_LL)", bad_nested == 0L, sprintf("%d violations", bad_nested))
ok("the ZC fit satisfies Cov(U,V) = 0", bad_cov == 0L, sprintf("%d violations", bad_cov))
ok("every fit is a probability vector", bad_sum == 0L,
   sprintf("max|sum-1| = %.1e, min(pi) = %.1e", worst_sum, worst_neg))
ok("the LL gradient is small", maxgrad < 1e-3, sprintf("max %.1e", maxgrad))

## The LL log-likelihood is concave: the same solution from any starting value
worst <- 0
for (b in 1:40) {
  r <- 4; u <- seq_len(r)
  y <- as.vector(rmultinom(1, 150, pop_pi(A_VECTORS$skew(r), 0.3, 0.5, u)))
  base <- fit_LL(y, r, u)$G2
  for (st in list(rep(0, r), c(2, -1, 1.5, 1), c(-2, 1, -1.5, -1)))
    worst <- max(worst, abs(fit_LL(y, r, u, start = st)$G2 - base))
}
ok("LL does not depend on the starting value", worst < 1e-6, sprintf("max %.1e", worst))

## With every pair in the first row the LL maximum is not attained: theta-hat
## diverges but the fit converges to the supremum, where G2(H_LL) is zero.
fb <- fit_LL(c(23, 13, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), 5, 1:5)
ok("LL on a boundary table: the supremum is returned",
   fb$G2 < 1e-6 && (fb$theta < 1e-6 || fb$theta > 1e6) && fb$gradient_max < 1e-6,
   sprintf("G2 %.1e, theta %.1e", fb$G2, fb$theta))

## (c) the public entry point ------------------------------------------
cat("\n-- (c) the entry point analyse_intraclass --\n")
ref <- analyse_intraclass(c(4, 17, 8, 21, 20, 6))$G2
M <- matrix(NA, 3, 3); M[1, ] <- c(4, 17, 8); M[2, 2:3] <- c(21, 20); M[3, 3] <- 6
ok("r inferred from the upper-triangular vector", ref[["I"]] > 0)
ok("upper-triangular matrix, lower triangle NA", isTRUE(all.equal(analyse_intraclass(M)$G2, ref)))
M0 <- M; M0[is.na(M0)] <- 0
ok("upper-triangular matrix, lower triangle 0", isTRUE(all.equal(analyse_intraclass(M0)$G2, ref)))
FULL <- matrix(c(4, 10, 5, 7, 21, 12, 3, 8, 6), 3, 3, byrow = TRUE)
ok("a full table is folded automatically", isTRUE(all.equal(analyse_intraclass(FULL)$G2, ref)))
ok("G2(H_I) does not depend on the scores",
   abs(analyse_intraclass(c(4, 17, 8, 21, 20, 6), scores = c(1, 2, 7))$G2[["I"]]
       - ref[["I"]]) < 1e-9)
bad_in <- function(expr) inherits(tryCatch(expr, error = function(e) e), "error")
ok("a length that is not r(r+1)/2 is an error", bad_in(analyse_intraclass(c(1, 2, 3, 4))))
ok("the wrong number of scores is an error",
   bad_in(analyse_intraclass(c(4, 17, 8, 21, 20, 6), scores = 1:4)))
ok("non-monotone scores are an error",
   bad_in(analyse_intraclass(c(4, 17, 8, 21, 20, 6), scores = c(3, 2, 1))))

## Are the populations what they claim to be? ---------------------------
cat("\n-- (d) the populations of the simulation --\n")
okpop <- TRUE
for (r in c(3, 4, 5)) {
  u <- seq_len(r); a <- A_VECTORS$skew(r)
  tz <- t_for_ZC(a, 0.6, u)
  okpop <- okpop &&
    is_I(pop_pi(a, 0, 0, u), r, u) &&                                   # Case 1
    is_ZC(pop_pi(a, tz, 0.6, u), r, u) && !is_LL(pop_pi(a, tz, 0.6, u), r, u) &&  # Case 2
    is_LL(pop_pi(a, 0.2, 0, u), r, u) && !is_ZC(pop_pi(a, 0.2, 0, u), r, u) &&    # Case 3
    !is_ZC(pop_pi(a, 0.2, 0.6, u), r, u) && !is_LL(pop_pi(a, 0.2, 0.6, u), r, u)  # Case 4
}
ok("the four cases are as intended at r = 3, 4, 5", okpop)

cat(sprintf("\n%s (%d failed)\n", if (fails == 0L) "all checks passed" else "some checks failed", fails))
if (fails > 0L) quit(status = 1)
