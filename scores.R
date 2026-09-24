# =====================================================================
#  Sensitivity of the decomposition to the choice of scores.
#
#  G^2(H_I) does not involve the scores and is therefore invariant;
#  the two components into which Theorem 2 splits it are not. This
#  script quantifies the split under several monotone score systems.
# =====================================================================

source("intraclass.R")

#' Midrank (ridit-type) scores from the estimated common marginal.
midrank_scores <- function(y, r) {
  n <- sum(y)
  p <- row_col_total(y, r) / (2 * n)
  cum <- cumsum(p)
  (c(0, cum[-r]) + cum) / 2
}
#' Normal (van der Waerden) scores from the same marginal.
normal_scores <- function(y, r) stats::qnorm(midrank_scores(y, r))

SETS <- function(y, r) {
  out <- list()
  out[["equally spaced"]] <- seq_len(r)
  if (r == 3) {
    out[["(1,2,4)"]] <- c(1, 2, 4)
    out[["(1,3,4)"]] <- c(1, 3, 4)
    out[["(1,2,10)"]] <- c(1, 2, 10)
  } else if (r == 4) {
    out[["(1,2,3,5)"]] <- c(1, 2, 3, 5)
    out[["(1,3,4,5)"]] <- c(1, 3, 4, 5)
    out[["(1,2,3,10)"]] <- c(1, 2, 3, 10)
  }
  out[["midranks"]] <- midrank_scores(y, r)
  out[["normal"]]   <- normal_scores(y, r)
  out
}

report <- function(name, y, r) {
  cat(sprintf("\n=== %s  (r = %d, n = %d) ===\n", name, r, sum(y)))
  cat(sprintf("%-16s %9s %9s %9s %9s %9s %9s\n",
              "scores", "G2(ZC)", "p", "G2(LL)", "p", "G2(I)", "theta"))
  for (nm in names(SETS(y, r))) {
    u <- SETS(y, r)[[nm]]
    d <- decompose(y, r, u)
    cat(sprintf("%-16s %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f\n",
                nm, d$G2["ZC"], d$p["ZC"], d$G2["LL"], d$p["LL"], d$G2["I"], d$theta))
  }
}

report("Ishii (1960), paired learning", c(4, 17, 8, 21, 20, 6), 3)
report("Data 1", c(14, 85, 194, 65, 81, 272, 136, 263, 318, 42), 4)
report("Data 2", c(45, 85, 74, 33, 65, 93, 69, 62, 52, 24), 4)
