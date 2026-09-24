# =====================================================================
#  Probability tables used in the simulation study.
#
#  All tables come from the one family
#
#      pi*_ij  propto  exp( a_i + a_j + t u_i u_j + delta 1{i = j} ),
#
#  with pi_ii = pi*_ii and pi_ij = 2 pi*_ij for i < j. The matrix Pi_* is
#  symmetric by construction, so the full table p_ij = p_ji = pi_ij / 2 is
#  symmetric and therefore satisfies marginal homogeneity; the intraclass
#  probabilities are all that the sampling scheme uses.
#
#  Inside the family,
#      t = 0, delta = 0   ->  H_I   (and hence both H_ZC and H_LL)
#      delta = 0, t != 0  ->  H_LL holds with theta = e^t, so by Theorem 1
#                             H_ZC must fail
#      delta != 0         ->  H_LL fails; H_ZC holds only for the one value
#                             of t that sets the covariance to zero
# =====================================================================

if (!exists("cov_UV")) source("intraclass.R")

#' Cell probabilities of the family.
pop_pi <- function(a, t, delta, u) {
  r <- length(a)
  W <- exp(outer(a, a, "+") + t * outer(u, u) + delta * diag(r))
  W <- W / sum(W)
  star_to_pi(W, r)
}

#' Does the LL association model hold at these probabilities?
#' LL holds iff the local log odds ratios of Pi_*, divided by the
#' corresponding product of score differences, are all equal.
is_LL <- function(pi_vec, r, u, tol = 1e-8) {
  M <- log(pi_star(pi_vec, r))
  v <- c()
  for (i in seq_len(r - 1)) for (j in seq_len(r - 1))
    v <- c(v, (M[i, j] + M[i + 1, j + 1] - M[i, j + 1] - M[i + 1, j]) /
               ((u[i + 1] - u[i]) * (u[j + 1] - u[j])))
  diff(range(v)) < tol
}

#' Does the zero covariance model hold?
is_ZC <- function(pi_vec, r, u, tol = 1e-8) abs(cov_UV(pi_vec, r, u)) < tol

#' Does the independence model hold?
is_I <- function(pi_vec, r, u, tol = 1e-8)
  is_LL(pi_vec, r, u, tol) && is_ZC(pi_vec, r, u, tol)

#' Value of t that makes the covariance vanish, for given a and delta.
t_for_ZC <- function(a, delta, u, interval = c(-3, 3)) {
  r <- length(a)
  f <- function(t) cov_UV(pop_pi(a, t, delta, u), r, u)
  if (f(interval[1]) * f(interval[2]) > 0) return(NA_real_)
  stats::uniroot(f, interval, tol = 1e-12)$root
}

#' Population divergence per observation from a model, i.e. the limit of
#' G^2 / n. Obtained by fitting the model to the expected counts.
pop_divergence <- function(pi_vec, r, u, model = c("I", "ZC", "LL")) {
  model <- match.arg(model)
  y <- pi_vec * 1e6
  f <- switch(model, I = fit_I(y, r), ZC = fit_ZC(y, r, u), LL = fit_LL(y, r, u))
  if (is.null(f)) return(NA_real_)
  f$G2 / 1e6
}

#' The marginal shapes used in the study.
#'
#' Two families are provided. `A_VECTORS` is the original one, in which the
#' log-weights are equally spaced over a fixed range; it is retained because
#' the size and power study of Section 5 was run with it.
#'
#' `a_entropy` is the graded family used for the study of the partition. A
#' label such as "0.85" there means the normalised Shannon entropy
#' H(p) / log r, which equals one at the uniform distribution and decreases as
#' the distribution concentrates. Fixing the entropy rather than the range
#' makes the marginal comparable across r: under `A_VECTORS` the same name
#' denotes entropy 0.942 at r = 3 but 0.970 at r = 5, so the number of
#' categories and the degree of concentration are confounded.
A_VECTORS <- list(
  flat   = function(r) rep(0, r),
  skew   = function(r) seq(0, -0.9, length.out = r),
  sparse = function(r) seq(0, -1.8, length.out = r)
)

#' Normalised Shannon entropy of the marginal implied by log-weights `a`.
norm_entropy <- function(a) {
  p <- exp(a - max(a)); p <- p / sum(p)
  -sum(p * log(p)) / log(length(p))
}

#' Log-weights whose marginal has the given normalised entropy.
#' @param target entropy in (0, 1]; 1 gives the uniform distribution.
a_entropy <- function(target, r) {
  if (target >= 1 - 1e-12) return(rep(0, r))
  f <- function(cc) norm_entropy(seq(0, -cc, length.out = r)) - target
  cc <- stats::uniroot(f, c(1e-8, 40), tol = 1e-10)$root
  seq(0, -cc, length.out = r)
}

#' The five levels used in the study of the partition.
#' New levels must be appended, never inserted: `make_grid` varies `shape`
#' slowest, so the position of a level fixes the condition numbers of its
#' block, and inserting one would silently reassign the files already written.
ENTROPY_LEVELS <- c(1.00, 0.95, 0.85, 0.75, 0.50)
