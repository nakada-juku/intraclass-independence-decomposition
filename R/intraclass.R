# =====================================================================
#  Decomposition of the independence model for two-way intraclass
#  contingency tables with ordered categories
#
#  Core routines: models H_I, H_ZC, H_LL and their likelihood ratio tests.
#
#  An intraclass table is stored as the upper triangle (including the
#  diagonal) of an r x r array, in the order
#      (1,1),(1,2),...,(1,r),(2,2),...,(2,r),...,(r,r),
#  so a table of order r is a vector of length K = r(r+1)/2.
# =====================================================================

## ---- indexing -------------------------------------------------------

#' Row/column indices of the upper triangle, in storage order.
upper_index <- function(r) {
  ij <- which(upper.tri(matrix(0, r, r), diag = TRUE), arr.ind = TRUE)
  ij[order(ij[, "row"], ij[, "col"]), , drop = FALSE]
}

#' Intraclass vector -> the symmetric matrix Pi_* of the paper.
#' Pi_*[i,i] = pi_ii and Pi_*[i,j] = pi_ij / 2 for i != j.
pi_star <- function(pi_vec, r) {
  idx <- upper_index(r)
  M <- matrix(0, r, r)
  for (k in seq_len(nrow(idx))) {
    i <- idx[k, 1]; j <- idx[k, 2]
    if (i == j) M[i, i] <- pi_vec[k] else { M[i, j] <- pi_vec[k] / 2; M[j, i] <- M[i, j] }
  }
  M
}

#' The symmetric matrix Pi_* -> intraclass vector.
star_to_pi <- function(M, r) {
  idx <- upper_index(r)
  apply(idx, 1, function(v) if (v[1] == v[2]) M[v[1], v[1]] else 2 * M[v[1], v[2]])
}

#' Fold a full r x r table into its intraclass version.
fold_table <- function(A) {
  r <- nrow(A)
  idx <- upper_index(r)
  apply(idx, 1, function(v) if (v[1] == v[2]) A[v[1], v[1]] else A[v[1], v[2]] + A[v[2], v[1]])
}

#' n_i^* = n_{i+} + n_{+i}, computed from the intraclass counts:
#'   n_i^* = 2 n_ii + sum_{j<i} n_ji + sum_{j>i} n_ij,
#' a diagonal pair contributing twice (Section 2 of the paper).
n_star <- function(y, r) {
  idx <- upper_index(r)
  tot <- numeric(r)
  for (k in seq_len(nrow(idx))) {
    i <- idx[k, 1]; j <- idx[k, 2]
    if (i == j) tot[i] <- tot[i] + 2 * y[k]
    else { tot[i] <- tot[i] + y[k]; tot[j] <- tot[j] + y[k] }
  }
  tot
}

## ---- the three models -----------------------------------------------

#' Covariance of (U, V) under the intraclass probabilities, assuming MH.
#'   Cov(U,V) = u' (Pi_* - Pi_* 1 1' Pi_*) u
cov_UV <- function(pi_vec, r, u) {
  M <- pi_star(pi_vec, r)
  as.numeric(t(u) %*% (M - M %*% rep(1, r) %*% t(rep(1, r)) %*% M) %*% u)
}

#' Likelihood ratio statistic 2 * sum n log(n / m); cells with n = 0 drop out.
G2_stat <- function(y, m) {
  ok <- y > 0
  2 * sum(y[ok] * log(y[ok] / m[ok]))
}

#' MLE under the independence model H_I (closed form).
#'   alpha_i = (n_{i+} + n_{+i}) / (2n)
fit_I <- function(y, r) {
  n <- sum(y)
  alpha <- n_star(y, r) / (2 * n)
  idx <- upper_index(r)
  pi_hat <- apply(idx, 1, function(v)
    if (v[1] == v[2]) alpha[v[1]]^2 else 2 * alpha[v[1]] * alpha[v[2]])
  m <- n * pi_hat
  list(pi = pi_hat, m = m, alpha = alpha,
       G2 = G2_stat(y, m), df = r * (r - 1) / 2)
}

#' Cell probabilities of the LL association model from (a, t),
#' where alpha_i = exp(a_i) and theta = exp(t). The normalisation
#' sum_{i,j} alpha_i alpha_j theta^{u_i u_j} = 1 is imposed by dividing
#' through, so (a, t) is unconstrained; a_1 is fixed at 0 for identifiability.
ll_pi <- function(par, r, u) {
  a <- c(0, par[seq_len(r - 1)])
  t <- par[r]
  W <- exp(outer(a, a, "+") + t * outer(u, u))   # alpha_i alpha_j theta^{u_i u_j}
  W <- W / sum(W)
  star_to_pi(W, r)
}

#' MLE under the LL association model H_LL.
#'
#' Writing pi*_ij proportional to exp(a_i + a_j + t u_i u_j), the model is an
#' exponential family and the log-likelihood
#'
#'   l(a, t) = sum_{i<=j} y_ij (a_i + a_j + t u_i u_j) + const
#'                - n log sum_{i,j} exp(a_i + a_j + t u_i u_j)
#'
#' is concave: the first term is linear and the log-sum-exp is convex. The
#' maximum is therefore unique once a_1 is fixed at zero, and one run of a
#' gradient method from any starting value finds it; no multi-start search is
#' needed. The score equations match the observed values of n_{i+} + n_{+i}
#' and of E[UV], as they must in an exponential family.
#'
#' @param start optional starting value for (a_2, ..., a_r, t).
fit_LL <- function(y, r, u, start = NULL) {
  n <- sum(y)
  # The scores are standardised internally. Both H_LL and H_ZC are invariant
  # under an increasing affine change u -> a u + b (Proposition 1), so the fit
  # is unaffected, but the exponent t u_i u_j is kept within range and large
  # scores no longer overflow. theta is not invariant and is converted back by
  # theta(u) = theta(u_std)^(1/s^2).
  s_u <- stats::sd(u); if (!is.finite(s_u) || s_u <= 0) s_u <- 1
  u <- (u - mean(u)) / s_u
  idx <- upper_index(r)
  tot <- n_star(y, r)          # n_{i+} + n_{+i}
  suv <- sum(y * u[idx[, 1]] * u[idx[, 2]])   # sum_{i<=j} y_ij u_i u_j
  UU  <- outer(u, u)

  pars <- function(par) list(a = c(0, par[seq_len(r - 1)]), t = par[r])
  star <- function(par) {             # the symmetric matrix Pi_*
    z <- pars(par)
    W <- exp(outer(z$a, z$a, "+") + z$t * UU)
    W / sum(W)
  }
  negll <- function(par) {
    z <- pars(par)
    E <- outer(z$a, z$a, "+") + z$t * UU
    m <- max(E)
    lS <- m + log(sum(exp(E - m)))    # log-sum-exp, computed stably
    -(sum(tot * z$a) / 2 * 2 - 0)     # placeholder, replaced below
  }
  # sum_{i<=j} y_ij (a_i + a_j) = sum_i a_i (n_{i+} + n_{+i}) / 1
  negll <- function(par) {
    z <- pars(par)
    E <- outer(z$a, z$a, "+") + z$t * UU
    m <- max(E); lS <- m + log(sum(exp(E - m)))
    -(sum(tot * z$a) + z$t * suv - n * lS)
  }
  neggr <- function(par) {
    P <- star(par)                    # Pi_*
    g_a <- tot - 2 * n * rowSums(P)   # d l / d a_k
    g_t <- suv - n * sum(P * UU)      # d l / d t
    -c(g_a[-1], g_t)
  }
  s0 <- if (!is.null(start)) start else {
    al <- pmax(fit_I(y, r)$alpha, 1e-12)
    c(log(al[-1]) - log(al[1]), 0)
  }
  o <- tryCatch(optim(s0, negll, neggr, method = "BFGS",
                      control = list(maxit = 2000, reltol = 1e-15)),
                error = function(e) NULL)
  if (is.null(o)) return(NULL)
  # The log-likelihood is concave, so the maximum is unique. The gradient left
  # by the BFGS stopping rule is driven to machine precision by damped Newton
  # steps with the analytic Hessian: log S is a partition function, so its
  # Hessian is the covariance of the sufficient statistic under Pi_*.
  hess <- function(par) {
    P <- star(par)
    rs <- rowSums(P)                       # sum_j pi*_kj
    Ea <- 2 * rs[-1]                       # E[s_k], k = 2..r
    Et <- sum(P * UU)                      # E[s_t]
    M  <- 2 * P[-1, -1, drop = FALSE]      # E[s_k s_l], k != l
    diag(M) <- 2 * rs[-1] + 2 * diag(P)[-1]
    Ct <- 2 * as.vector((P * UU)[-1, , drop = FALSE] %*% rep(1, r))  # E[s_k s_t]
    Ett <- sum(P * UU^2)
    H <- rbind(cbind(M - outer(Ea, Ea), Ct - Ea * Et),
               c(Ct - Ea * Et, Ett - Et^2))
    n * H
  }
  par <- o$par; val <- negll(par)
  for (it in seq_len(30)) {
    g <- neggr(par)
    if (max(abs(g)) < 1e-10) break
    st <- tryCatch(solve(hess(par) + diag(1e-10, r), -g), error = function(e) NULL)
    if (is.null(st)) break
    tt <- 1; moved <- FALSE
    repeat {
      cand <- par + tt * st; v <- negll(cand)
      if (is.finite(v) && v <= val) { par <- cand; val <- v; moved <- TRUE; break }
      tt <- tt / 2; if (tt < 1e-12) break
    }
    if (!moved) break
  }
  o$par <- par
  pi_hat <- star_to_pi(star(o$par), r)
  m <- n * pi_hat
  list(pi = pi_hat, m = m, theta = exp(o$par[r])^(1 / s_u^2),
       alpha = exp(c(0, o$par[seq_len(r - 1)])),
       G2 = G2_stat(y, m), df = (r - 2) * (r + 1) / 2,
       converged = o$convergence == 0,
       gradient_max = max(abs(neggr(o$par))))
}

#' Coefficients that express the covariance in the intraclass probabilities.
#' Writing uu_k = u_i u_j and ubar_k = (u_i + u_j)/2 for the kth cell (i,j),
#'   E[UV] = uu' pi,   mu = E[U] = ubar' pi,   Cov(U,V) = uu' pi - (ubar' pi)^2.
#' The symbol g is deliberately not used here: in the proof of Theorem 2 the
#' paper writes g_ij = u_i u_j - (u_i + u_j) mu, which is a different quantity.
cov_coef <- function(r, u) {
  idx <- upper_index(r)
  list(uu = u[idx[, 1]] * u[idx[, 2]],
       ubar = (u[idx[, 1]] + u[idx[, 2]]) / 2)
}

#' MLE under the zero covariance model H_ZC (numerical).
#'
#' The constraint Cov(U,V) = g'pi - (s'pi)^2 = 0 is quadratic and not convex,
#' so a Lagrangian on it can have several stationary points and a penalty
#' formulation can jump between branches. It becomes linear once the value of
#' the mean m = s'pi is held fixed:
#'
#'     Cov = 0   <=>   there is an m with   s'pi = m   and   g'pi = m^2 .
#'
#' For fixed m the problem is therefore the maximisation of a concave
#' log-likelihood over the simplex subject to two further linear constraints,
#' which has a unique solution; it is solved here in its three-dimensional
#' convex dual. The profile log-likelihood is then maximised over the scalar m.
#'
#' A cell with no observations may still receive positive probability under
#' this model: it contributes nothing to the likelihood but does relax the two
#' moment constraints, so the maximum is not in general attained at zero there.
#' This is handled by the barrier described inside the function.
fit_ZC <- function(y, r, u, tol = 1e-9, ngrid = 60) {
  n <- sum(y)
  # Standardised for the same reason as in fit_LL. The constraint Cov(U,V) = 0
  # is invariant under an affine change, so the fit is unaffected; mu is
  # converted back to the original scale.
  s_u <- stats::sd(u); if (!is.finite(s_u) || s_u <= 0) s_u <- 1
  m_u <- mean(u); u <- (u - m_u) / s_u
  p0 <- y / n
  h0 <- cov_UV(p0, r, u)
  if (abs(h0) < tol * max(1, abs(h0)))
    return(list(pi = p0, m = y, G2 = 0, df = 1, converged = TRUE,
                mu_hat = s_u * sum(cov_coef(r, u)$ubar * p0) + m_u, cov_at_opt = h0))

  cc <- cov_coef(r, u); uu <- cc$uu; ubar <- cc$ubar
  K <- length(y); sup <- which(y > 0); emp <- which(y == 0)
  A <- rbind(rep(1, K), ubar, uu)                   # 3 x K, all cells

  # For fixed m: maximise sum_{y_k>0} y_k log pi_k subject to A pi = (1,m,m^2)'
  # and pi >= 0. Cells with y_k = 0 do not enter the objective but may still
  # carry mass, because they relax the two moment constraints; the Karush--Kuhn
  # --Tucker conditions give pi_k = y_k / (lambda' A_k) on the support and
  # lambda' A_k >= 0 off it. The dual is therefore minimised subject to those
  # inequalities, which is done here by a logarithmic barrier: with weights
  # w_k = y_k on the support and w_k = mu off it, the stationary point of
  #     -sum_k w_k log(lambda' A_k) + lambda' b
  # satisfies A pi = b exactly for pi_k = w_k / (lambda' A_k), and letting
  # mu decrease drives the off-support masses to their optimal values.
  # Each barrier subproblem is convex and is solved by a damped Newton step
  # with Armijo backtracking that keeps lambda' A_k > 0.
  # Barrier parameter. Named apart from mu, which is E[U] in the paper.
  BARRIER <- c(1, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6)
  inner <- function(mu, lam0 = NULL, tol_grad = 1e-8, maxit = 40) {
    b <- c(1, mu, mu^2)
    lam <- if (is.null(lam0)) c(n, 0, 0) else lam0   # c(n,0,0) gives pi = y/n
    if (any(as.vector(lam %*% A) <= 0)) lam <- c(n, 0, 0)
    bars <- if (length(emp)) BARRIER else 0
    for (bar in bars) {
      w <- y; w[emp] <- bar
      Dfun <- function(lm) {
        v <- as.vector(lm %*% A)
        if (any(v <= 0)) return(Inf)
        -sum(w * log(v)) + sum(lm * b)
      }
      Dc <- Dfun(lam); if (!is.finite(Dc)) return(NULL)
      # The gradient is b - A pi, so its smallness is exactly feasibility.
      # Stalling of the objective is not a stopping rule: it would leave the
      # constraints only loosely satisfied.
      for (it in seq_len(maxit)) {
        v <- as.vector(lam %*% A)
        grad <- as.vector(-(A %*% (w / v))) + b
        if (max(abs(grad)) < tol_grad) break
        H <- (A * rep(w / v^2, each = 3)) %*% t(A)
        step <- tryCatch(solve(H + diag(1e-14, 3), -grad), error = function(e) NULL)
        if (is.null(step)) return(NULL)
        slope <- sum(grad * step)          # negative for a descent direction
        tt <- 1
        repeat {                            # Armijo backtracking inside the cone
          lam_new <- lam + tt * step
          Dn <- Dfun(lam_new)
          if (is.finite(Dn) && Dn <= Dc + 1e-4 * tt * slope) break
          tt <- tt / 2
          if (tt < 1e-14) return(NULL)
        }
        stalled <- abs(Dc - Dn) < 1e-16 * max(1, abs(Dc))
        lam <- lam_new; Dc <- Dn
        if (stalled && tt < 1e-8) break   # give up only when no step helps
      }
    }
    v <- as.vector(lam %*% A)
    if (any(v <= 0)) return(NULL)
    w <- y; if (length(emp)) w[emp] <- BARRIER[length(BARRIER)]
    p <- w / v
    if (abs(sum(p) - 1) > 1e-6 ||
        abs(sum(ubar * p) - mu) > 1e-6 ||
        abs(sum(uu * p) - mu^2) > 1e-6) return(NULL)
    attr(p, "lam") <- lam
    p
  }
  NEG <- -1e300
  loglik <- function(p) if (is.null(p)) NEG else sum(y[sup] * log(p[sup]))

  # Search over mu, which lies in the range of ubar. The profile
  # log-likelihood need not be unimodal, so a grid is scanned first and the
  # best two local maxima are then refined.
  lo <- min(ubar); hi <- max(ubar)  # empty cells may carry mass: use all cells
  if (lo >= hi) return(NULL)
  fI <- fit_I(y, r)
  mus <- sort(unique(c(seq(lo + 1e-6 * (hi - lo), hi - 1e-6 * (hi - lo),
                          length.out = ngrid),
                      sum(ubar * fI$pi))))            # mu of the independence fit
  ng <- length(mus)
  # sweep outwards from the grid point nearest the unrestricted mean, carrying
  # the dual solution forward as the starting value for the next m
  mu0 <- sum(ubar * p0); j0 <- which.min(abs(mus - mu0))
  ps <- vector("list", ng)
  lam <- NULL
  for (j in j0:ng) { v <- inner(mus[j], lam); ps[j] <- list(v)
                     if (!is.null(v)) lam <- attr(v, "lam") }
  lam <- if (is.null(ps[[j0]])) NULL else attr(ps[[j0]], "lam")
  if (j0 > 1) for (j in (j0 - 1):1) { v <- inner(mus[j], lam); ps[j] <- list(v)
                                      if (!is.null(v)) lam <- attr(v, "lam") }
  lls <- vapply(ps, loglik, numeric(1))
  if (all(lls <= NEG)) {
    return(list(pi = fI$pi, m = n * fI$pi, G2 = fI$G2, df = 1,
                converged = FALSE, mu_hat = s_u * sum(ubar * fI$pi) + m_u,
                cov_at_opt = cov_UV(fI$pi, r, u)))
  }
  k <- which.max(lls)
  best_mu <- mus[k]; best_p <- ps[[k]]; best_ll <- lls[k]
  # refine around the best grid points and around any other local maximum
  top <- order(lls, decreasing = TRUE)[seq_len(min(3, sum(lls > NEG)))]
  peaks <- which(lls > NEG &
                 c(TRUE, diff(lls) > 0) & c(diff(lls) < 0, TRUE))
  for (j in unique(c(k, top, peaks))) {
    a1 <- mus[max(1, j - 2)]; a2 <- mus[min(ng, j + 2)]
    if (a2 <= a1) next
    ref <- tryCatch(optimize(function(mu) loglik(inner(mu)), c(a1, a2),
                             maximum = TRUE, tol = 1e-11),
                    error = function(e) NULL)
    if (!is.null(ref) && ref$objective > best_ll) {
      cand <- inner(ref$maximum)
      if (!is.null(cand)) { best_mu <- ref$maximum; best_p <- cand; best_ll <- ref$objective }
    }
  }

  # The independence fit satisfies Cov(U,V) = 0 exactly, because Pi_* = alpha
  # alpha' with 1'alpha = 1, so it is feasible here; keeping it as a candidate
  # guarantees that the fit is never worse than under the nested model H_I.
  # The grid is scanned at a loose tolerance and only the selected mu is
  # re-solved tightly.
  polished <- inner(best_mu, attr(best_p, "lam"), tol_grad = 1e-13, maxit = 200)
  if (!is.null(polished) && loglik(polished) >= best_ll - 1e-9) best_p <- polished

  if (loglik(fI$pi) > best_ll) { best_p <- fI$pi; best_mu <- sum(ubar * fI$pi) }

  m <- n * best_p
  list(pi = best_p, m = m, G2 = G2_stat(y, m), df = 1,
       converged = TRUE, mu_hat = s_u * best_mu + m_u,
       cov_at_opt = cov_UV(best_p, r, u))
}

## ---- the three tests together ---------------------------------------

#' Fit H_I, H_ZC and H_LL to one intraclass table and return the
#' likelihood ratio statistics, their degrees of freedom and p-values,
#' together with the discrepancy of the partition of Theorem 2.
#'
#' @param y intraclass counts, upper triangle in storage order
#' @param r number of categories
#' @param u scores; defaults to 1, ..., r
decompose <- function(y, r, u = seq_len(r)) {
  fI <- fit_I(y, r); fZ <- fit_ZC(y, r, u); fL <- fit_LL(y, r, u)
  if (is.null(fZ) || is.null(fL)) return(NULL)
  # G2 and df are the G^2(H_M) and degrees of freedom of the paper, named for
  # H_I, H_ZC and H_LL.
  G2 <- c(I = fI$G2, ZC = fZ$G2, LL = fL$G2)
  df <- c(I = fI$df, ZC = fZ$df, LL = fL$df)
  # H_I is nested in both H_ZC and H_LL, so G2(H_I) must dominate the other
  # two up to rounding. A violation means one of the numerical fits failed.
  tolv <- 1e-6 * max(1, G2["I"])
  valid <- unname(G2["ZC"] <= G2["I"] + tolv && G2["LL"] <= G2["I"] + tolv)
  list(G2 = G2, df = df,
       p = stats::pchisq(G2, df, lower.tail = FALSE),
       theta = fL$theta,
       # D is the D of the paper: G^2(H_I) - {G^2(H_ZC) + G^2(H_LL)}
       D = unname(G2["I"] - G2["ZC"] - G2["LL"]),
       valid = valid,
       converged = isTRUE(fZ$converged) && isTRUE(fL$converged),
       fits = list(I = fI, ZC = fZ, LL = fL))
}

## ---- public entry point ---------------------------------------------

#' Coerce user input to the intraclass storage vector, inferring r.
#'
#' @param y one of
#'   (a) a numeric vector holding the upper triangle, including the
#'       diagonal, in the order (1,1),(1,2),...,(1,r),(2,2),...,(r,r);
#'   (b) an r x r matrix whose upper triangle holds the intraclass counts,
#'       the lower triangle being zero or NA;
#'   (c) a full r x r table of ordered pairs, which is folded.
#' @param fold `TRUE` to fold, `FALSE` to read the upper triangle as it
#'   stands, or `NULL` (default) to decide from the lower triangle.
as_intraclass <- function(y, fold = NULL) {
  if (is.matrix(y) || is.data.frame(y)) {
    A <- as.matrix(y)
    if (nrow(A) != ncol(A)) stop("the matrix must be square")
    r <- nrow(A)
    lower <- A[lower.tri(A)]
    if (is.null(fold)) fold <- any(!is.na(lower) & lower != 0)
    if (fold) {
      A[is.na(A)] <- 0
      return(list(y = fold_table(A), r = r, folded = TRUE))
    }
    A[is.na(A)] <- 0
    idx <- upper_index(r)
    return(list(y = A[idx], r = r, folded = FALSE))
  }
  y <- as.numeric(y)
  K <- length(y)
  r <- (-1 + sqrt(1 + 8 * K)) / 2
  if (abs(r - round(r)) > 1e-8)
    stop("length ", K, " is not of the form r(r+1)/2")
  list(y = y, r = round(r), folded = FALSE)
}

#' Decomposition of the independence model for one intraclass table.
#'
#' Give it a table and it returns everything: the likelihood ratio
#' statistics for H_I, H_ZC and H_LL with their degrees of freedom and
#' p-values, the estimate of theta, the fitted expected frequencies, the
#' discrepancy of the partition of Theorem 2, and the combined p-values.
#'
#' @param y      the table; see `as_intraclass`.
#' @param scores scores u_1 < ... < u_r. Defaults to 1, ..., r. The number
#'   of categories is inferred from `y` and need not be supplied.
#' @param fold   see `as_intraclass`.
#'
#' @examples
#'   ## Ishii (1960), progress of paired learning, r = 3, n = 76
#'   analyse_intraclass(c(4, 17, 8, 21, 20, 6))
analyse_intraclass <- function(y, scores = NULL, fold = NULL) {
  z <- as_intraclass(y, fold)
  r <- z$r; yy <- z$y
  u <- if (is.null(scores)) seq_len(r) else as.numeric(scores)
  if (length(u) != r) stop("expected ", r, " scores")
  if (any(diff(u) <= 0)) stop("the scores must be strictly increasing")
  d <- decompose(yy, r, u)
  if (is.null(d)) stop("the fit failed")
  p <- d$p
  comb <- c(Bonferroni = min(1, 2 * min(p["ZC"], p["LL"])),
            Simes      = { q <- sort(c(p["ZC"], p["LL"])); unname(min(2 * q[1], q[2])) },
            Fisher     = stats::pchisq(-2 * sum(log(pmax(c(p["ZC"], p["LL"]), 1e-300))),
                                       df = 4, lower.tail = FALSE))
  structure(list(r = r, n = sum(yy), scores = u, counts = yy, folded = z$folded,
                 G2 = d$G2, df = d$df, p = p, theta = d$theta,
                 partition_gap = d$D, combined_p = comb,
                 expected = lapply(d$fits, function(f) f$m),
                 valid = d$valid, converged = d$converged),
            class = "intraclass_decomposition")
}

#' @export
print.intraclass_decomposition <- function(x, digits = 3, ...) {
  cat(sprintf("Intraclass contingency table:  r = %d,  n = %g%s\n",
              x$r, x$n, if (x$folded) "  (folded from a full table)" else ""))
  cat("Scores:", paste(format(x$scores, digits = digits), collapse = ", "), "\n\n")
  tab <- data.frame(G2 = x$G2, df = x$df, p.value = x$p)
  rownames(tab) <- c("H_I  (independence)", "H_ZC (zero covariance)",
                     "H_LL (linear-by-linear)")
  print(round(tab, digits))
  cat(sprintf("\ntheta.hat = %s\n", format(x$theta, digits = digits)))
  cat(sprintf("G2(H_I) - {G2(H_ZC) + G2(H_LL)} = %s\n",
              format(x$partition_gap, digits = 2)))
  cat("\nCombined p-values for H_I from the two components:\n  ")
  cat(paste(sprintf("%s %s", names(x$combined_p),
                    format(x$combined_p, digits = digits)), collapse = "   "), "\n")
  if (!x$valid || !x$converged)
    cat("\nNote: a fit did not converge, or the fits are not consistent\n")
  invisible(x)
}
