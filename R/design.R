# =====================================================================
#  The design of the simulation study, and the code that runs one
#  condition of it. Sourced by run_simulation.R.
# =====================================================================

ALPHA <- 0.05

## ---- the alternatives ------------------------------------------------
## Each row names a population: `case` labels which of H_ZC and H_LL holds,
## `delta` is the diagonal parameter and `t` the log of theta; for Case 2
## the value of t is solved for so that the covariance vanishes.
SETTINGS <- data.frame(
  case  = c("1", "2a", "2b", "3a", "3b", "4a", "4b"),
  delta = c(  0,  0.3,  0.6,    0,    0,  0.4,  0.8),
  t     = c(  0,   NA,   NA,  0.1,  0.2, 0.15, 0.30),
  stringsAsFactors = FALSE
)

make_pop <- function(setting, r, shape, u) {
  a <- if (grepl("^H[0-9.]+$", shape)) a_entropy(as.numeric(sub("^H", "", shape)), r)
       else A_VECTORS[[shape]](r)
  tt <- if (is.na(setting$t)) t_for_ZC(a, setting$delta, u) else setting$t
  if (is.na(tt)) return(NULL)
  list(pi = pop_pi(a, tt, setting$delta, u), t = tt, a = a)
}

## ---- combining the two component p-values ----------------------------
p_bonferroni <- function(p1, p2) min(1, 2 * min(p1, p2))
p_simes      <- function(p1, p2) { p <- sort(c(p1, p2)); min(2 * p[1], p[2]) }
p_fisher     <- function(p1, p2)
  stats::pchisq(-2 * (log(max(p1, 1e-300)) + log(max(p2, 1e-300))),
                df = 4, lower.tail = FALSE)

## ---- one condition ----------------------------------------------------
run_condition <- function(i, GRID, B, u_fun = seq_len, seed_base = 0L) {
  g <- GRID[i, ]
  r <- g$r; n <- g$n; u <- u_fun(r)
  st <- SETTINGS[SETTINGS$case == g$case, ]
  pop <- make_pop(st, r, g$shape, u)
  if (is.null(pop)) return(NULL)
  pi0 <- pop$pi
  K <- r * (r + 1) / 2

  res <- matrix(NA_real_, B, 10)
  colnames(res) <- c("G2I", "G2ZC", "G2LL", "pI", "pZC", "pLL", "D", "zeros",
                     "valid", "conv")
  fails <- 0L
  for (b in seq_len(B)) {
    set.seed(seed_base + 10000L * i + b)
    y <- as.vector(rmultinom(1, n, pi0))
    d <- decompose(y, r, u)
    if (is.null(d)) { fails <- fails + 1L; next }
    res[b, ] <- c(d$G2["I"], d$G2["ZC"], d$G2["LL"],
                  d$p["I"], d$p["ZC"], d$p["LL"], d$D, sum(y == 0),
                  as.numeric(d$valid), as.numeric(d$converged))
  }
  ok <- !is.na(res[, "G2I"])
  R <- res[ok, , drop = FALSE]
  m <- nrow(R)
  if (m == 0) return(NULL)

  comb <- t(apply(R[, c("pZC", "pLL"), drop = FALSE], 1, function(p)
    c(bonf = p_bonferroni(p[1], p[2]),
      simes = p_simes(p[1], p[2]),
      fisher = p_fisher(p[1], p[2]))))

  data.frame(
    case = g$case, r = r, n = n, shape = g$shape,
    t_used = pop$t, delta = st$delta,
    d_I  = pop_divergence(pi0, r, u, "I"),
    d_ZC = pop_divergence(pi0, r, u, "ZC"),
    d_LL = pop_divergence(pi0, r, u, "LL"),
    B_done = m, fail_rate = 100 * fails / B,
    invalid_rate = 100 * mean(R[, "valid"] == 0),
    nonconv_rate = 100 * mean(R[, "conv"] == 0),
    # rejection rates (%)
    rej_LRT   = 100 * mean(R[, "pI"]  <= ALPHA),
    rej_ZC    = 100 * mean(R[, "pZC"] <= ALPHA),
    rej_LL    = 100 * mean(R[, "pLL"] <= ALPHA),
    rej_bonf  = 100 * mean(comb[, "bonf"]   <= ALPHA),
    rej_simes = 100 * mean(comb[, "simes"]  <= ALPHA),
    rej_fisher= 100 * mean(comb[, "fisher"] <= ALPHA),
    # accuracy of the partition of Theorem 2
    D_mean = mean(R[, "D"]), D_sd = sd(R[, "D"]),
    D_absmax = max(abs(R[, "D"])),
    D_q95 = unname(quantile(abs(R[, "D"]), 0.95)),
    # Signed quantiles, for plotting D itself on the vertical axis
    D_p05 = unname(quantile(R[, "D"], 0.05)),
    D_p25 = unname(quantile(R[, "D"], 0.25)),
    D_p50 = unname(quantile(R[, "D"], 0.50)),
    D_p75 = unname(quantile(R[, "D"], 0.75)),
    D_p95 = unname(quantile(R[, "D"], 0.95)),
    rel_D_med = median(abs(R[, "D"]) / pmax(R[, "G2I"], 1e-8)),
    # how well the statistic itself matches its reference distribution
    mean_G2I = mean(R[, "G2I"]), df_I = r * (r - 1) / 2,
    mean_G2ZC = mean(R[, "G2ZC"]), mean_G2LL = mean(R[, "G2LL"]),
    df_LL = (r - 2) * (r + 1) / 2,
    # sparseness
    zeros_mean = mean(R[, "zeros"]), K = K,
    any_zero = 100 * mean(R[, "zeros"] > 0),
    mcse = 100 * sqrt(0.05 * 0.95 / m),
    stringsAsFactors = FALSE
  )
}

## ---- the two designs --------------------------------------------------
#' `main`  : the four cases at four sample sizes, plus a sparse block.
#' `curve` : size and power as functions of n on a fine grid.
#' `order` : the rate at which the discrepancy D of Theorem 2 vanishes under
#'           H_I. Only Case 1 is needed, on a geometric grid of n wide enough
#'           for the slope of log s.d.(D) against log n to be read off.
make_grid <- function(design = c("main", "curve", "order", "order2", "order3")) {
  design <- match.arg(design)
  if (design == "main") {
    g <- expand.grid(case = SETTINGS$case, r = c(3, 4, 5),
                     n = c(50, 100, 200, 500), shape = "skew",
                     stringsAsFactors = FALSE)
    g <- rbind(g, expand.grid(case = c("1", "3a", "4a"), r = c(4, 5),
                              n = c(50, 100), shape = "sparse",
                              stringsAsFactors = FALSE))
  } else if (design == "order") {
    NS <- c(40, 70, 120, 200, 350, 600, 1000)
    g <- expand.grid(case = "1", r = c(3, 4, 5), n = NS,
                     shape = sprintf("H%.2f", ENTROPY_LEVELS),
                     stringsAsFactors = FALSE)
  } else if (design == "order3") {
    ## Carried far enough for e_min to exceed five at every marginal, which it
    ## does not under `order` and `order2` alone.
    NS <- c(10000)
    g <- expand.grid(case = "1", r = c(3, 4, 5), n = NS,
                     shape = sprintf("H%.2f", ENTROPY_LEVELS),
                     stringsAsFactors = FALSE)
  } else if (design == "order2") {
    ## Sample sizes added to those of `order`. Inserting an n into an existing
    ## design would shift the condition numbers, because n varies inside `shape`,
    ## and set.seed(10000 * i + b) would then point at a different stream.
    NS <- c(90, 150, 280, 2000, 4000)
    g <- expand.grid(case = "1", r = c(3, 4, 5), n = NS,
                     shape = sprintf("H%.2f", ENTROPY_LEVELS),
                     stringsAsFactors = FALSE)
  } else {
    NS <- c(seq(20, 200, by = 10), 250, 300, 400, 500)
    g <- expand.grid(case = c("1", "2b", "3b", "4a"), r = c(3, 4, 5),
                     n = NS, shape = "skew", stringsAsFactors = FALSE)
  }
  g
}
