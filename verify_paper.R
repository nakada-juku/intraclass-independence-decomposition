# =====================================================================
#  Check the manuscript against the stored results.
#
#      Rscript verify_paper.R [path/to/Maniscript_R1.tex]
#
#  Two kinds of check:
#    (a) every table in the manuscript is reproduced, row for row, by
#        make_tables.R;
#    (b) every number quoted in the prose of Section 5 is recomputed
#        from the results and compared with what is printed.
#
#  The point of (b) is that prose numbers are the ones that rot: a table
#  is regenerated when the data change, a sentence is not.
# =====================================================================

source("R/intraclass.R")
source("R/populations.R", chdir = TRUE)

args <- commandArgs(trailingOnly = TRUE)
TEX  <- if (length(args) >= 1) args[1] else "../Maniscript_R1.tex"

main  <- readRDS("results/main.rds")
order <- do.call(rbind, lapply(
  Filter(file.exists, sprintf("results/%s.rds", c("order", "order2", "order3"))), readRDS))
order$H    <- as.numeric(sub("^H", "", order$shape))
order$emin <- mapply(function(H, r, n)
  n * min(pop_pi(a_entropy(H, r), 0, 0, seq_len(r))), order$H, order$r, order$n)

nchk <- 0L; nbad <- 0L
ok <- function(label, got, want, tol = 0.0051) {
  nchk <<- nchk + 1L
  good <- !is.na(got) && abs(got - want) <= tol
  if (!good) nbad <<- nbad + 1L
  cat(sprintf("  %-50s printed %-10s computed %-10s %s\n", label, format(want),
              if (is.na(got)) "-" else sprintf("%.4g", got), if (good) "ok" else "MISMATCH"))
}

## ---- (a) the tables --------------------------------------------------
cat("-- tables reproduced by make_tables.R --\n")
tmp <- tempfile(fileext = ".tex")
system2("Rscript", c("make_tables.R", tmp), stdout = NULL, stderr = NULL)
datarow <- function(l) grepl("\\\\\\\\\\s*$", l) &&
  !grepl("rule|multicolumn|endhead|endfoot|endfirsthead|caption|Case &|data &|target &", l)
pull <- function(path, lbl, endpat) {
  s <- readLines(path, warn = FALSE)
  i <- grep(paste0("label\\{", lbl, "\\}"), s, fixed = FALSE)[1]
  j <- grep(endpat, s); j <- j[j > i][1]
  r <- s[(i + 1):(j - 1)]
  trimws(r[vapply(r, datarow, logical(1))])
}
for (lbl in c("Table_marginals", "Table_size", "Table_power", "Table_power_full", "Table_scores")) {
  endpat <- if (lbl == "Table_power_full") "end\\{longtable\\}" else "bottomrule"
  a <- pull(tmp, lbl, endpat); b <- pull(TEX, lbl, endpat)
  nchk <- nchk + 1L
  if (identical(a, b)) cat(sprintf("  %-50s %3d rows                   ok\n", lbl, length(a)))
  else { nbad <- nbad + 1L
    cat(sprintf("  %-50s generated %d / manuscript %d   MISMATCH\n", lbl, length(a), length(b))) }
}

## ---- (b) the prose ---------------------------------------------------
cat("\n-- numbers quoted in the prose of Section 5 --\n")
s1 <- main[main$case == "1", ]
ok("size of the ZC test, lower",  min(s1$rej_ZC), 4.7, 0.051)
ok("size of the ZC test, upper",  max(s1$rej_ZC), 5.9, 0.051)
ok("mean G2(H_I), r=5, n=50",     s1$mean_G2I[s1$r == 5 & s1$n == 50 & s1$shape == "skew"], 11.22)
c2 <- main[main$case %in% c("2a", "2b") & main$shape == "skew", ]
ok("Case 2, ZC size, lower",      min(c2$rej_ZC), 2.8, 0.051)
ok("Case 2, ZC size, upper",      max(c2$rej_ZC), 6.3, 0.051)
ok("largest mean of D over alternatives", max(main$D_mean), 15.1, 0.051)

sl <- c(); r2 <- c()
for (sh in unique(order$shape)) for (rr in sort(unique(order$r))) {
  d <- order[order$shape == sh & order$r == rr, ]
  f <- stats::lm(log(D_sd) ~ log(n), d)
  sl <- c(sl, stats::coef(f)[2]); r2 <- c(r2, summary(f)$r.squared)
}
ok("smallest fitted exponent", min(sl), SLOPE_LO <- round(min(sl), 2), 0.006)
ok("largest fitted exponent",  max(sl), SLOPE_HI <- round(max(sl), 2), 0.006)
cat(sprintf("  %-50s %s\n", "(exponents, for the text)",
            sprintf("%.3f to %.3f, mean %.3f, min R2 %.3f", min(sl), max(sl), mean(sl), min(r2))))

for (m in c(1, 2, 5, 10)) {
  z <- order[order$emin >= m, ]
  cat(sprintf("  %-50s %3d conditions, largest 95%% point %.3f, %d marginals\n",
              sprintf("e_min >= %d", m), nrow(z), max(z$D_q95), length(unique(z$shape))))
}
mx <- max(order$D_q95[order$emin >= 5])
w  <- order[order$emin >= 5, ][which.max(order$D_q95[order$emin >= 5]), ]
cat(sprintf("  %-50s %s r=%d n=%d\n", "worst condition with e_min >= 5", w$shape, w$r, w$n))
cat(sprintf("  %-50s %.1f%% at r=3, %.1f%% at r=5\n", "that bound as a share of E G2(H_I)",
            mx / 3 * 100, mx / 10 * 100))
cat(sprintf("  %-50s e_min %.3f / n %.3f / n per cell %.3f\n", "R2 of log 95%% point on",
            summary(stats::lm(log(D_q95) ~ log(emin), order))$r.squared,
            summary(stats::lm(log(D_q95) ~ log(n), order))$r.squared,
            summary(stats::lm(log(D_q95) ~ log(n / (r * (r + 1) / 2)), order))$r.squared))
cat(sprintf("  %-50s %.2f%% overall, %.2f%% where e_min >= 1\n", "non-convergence of the ZC fit",
            max(order$nonconv_rate), max(order$nonconv_rate[order$emin >= 1])))

cat(sprintf("\n%d checks, %d mismatches\n", nchk, nbad))
if (nbad > 0L) quit(status = 1L)
