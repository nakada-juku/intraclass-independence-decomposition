# =====================================================================
#  Every table of the paper, generated from the stored results.
#
#  Run from this directory, after collect_results.R:
#      Rscript make_tables.R [outfile]
#  The default output is tables.tex, which holds all five tables in the
#  order in which they appear in the paper.
# =====================================================================

source("R/intraclass.R")
source("R/populations.R", chdir = TRUE)
source("R/design.R")

args <- commandArgs(trailingOnly = TRUE)
OUT  <- if (length(args) >= 1) args[1] else "tables.tex"

main  <- readRDS("results/main.rds")
order <- do.call(rbind, lapply(
  Filter(file.exists, sprintf("results/%s.rds", c("order", "order2", "order3"))), readRDS))

con <- file(OUT, "w"); w <- function(...) writeLines(paste0(...), con)
LEVELS <- sprintf("%.2f", ENTROPY_LEVELS)

## ---- Table: the graded family of marginals --------------------------
w("%%%% TABLE: marginals")
w("\\begin{table}[htbp]")
w("\\centering")
w("\\caption{The graded family of marginal distributions of Subsection~\\ref{subsec_design},",
  " with the normalised entropy \\eqref{Eq_Entropy} attained in each case}")
w("\\label{Table_marginals}")
w("\\begin{tabular*}{\\textwidth}{@{\\extracolsep\\fill}llr} \\toprule")
w("target & $\\bm{p}$ & $\\mathcal{H}(\\bm{p})/\\log r$ \\\\ \\midrule")
for (r in c(3, 4, 5)) {
  w(sprintf("\\multicolumn{3}{l}{\\textit{$r=%d$}}\\\\", r))
  for (tg in ENTROPY_LEVELS) {
    a <- a_entropy(tg, r); p <- exp(a) / sum(exp(a))
    w(sprintf("$%.2f$ & $(%s)$ & $%.4f$ \\\\", tg,
              paste(sprintf("%.3f", p), collapse = ",\\,"), norm_entropy(a)))
  }
  if (r < 5) w("\\midrule")
}
w("\\bottomrule\\end{tabular*}\\end{table}")

## ---- Table: size ----------------------------------------------------
w("")
w("%%%% TABLE: size")
w("\\begin{table}[htbp]")
w("\\centering")
w("\\caption{Empirical size (\\%) of the tests at the nominal $5\\%$ level under $H_{I}$,",
  " from $5{,}000$ replications per condition, with the percentage of samples containing",
  " at least one empty cell in the last column}")
w("\\label{Table_size}")
w("\\begin{tabular*}{\\textwidth}{@{\\extracolsep\\fill}llrrrrrrr} \\toprule")
w("$\\bm{p}$ & $r$ & $n$ & $G^{2}(H_{I})$ & $G^{2}(H_{ZC})$ & $G^{2}(H_{LL})$",
  " & Bonferroni & Simes & Fisher \\\\ \\midrule")
d <- main[main$case == "1", ]; d <- d[order(-xtfrm(d$shape), d$r, d$n), ]; prev <- ""
for (i in seq_len(nrow(d))) {
  g <- d[i, ]
  sh <- if (g$shape == prev) "" else if (g$shape == "skew") "moderate" else "sparse"
  prev <- g$shape
  w(sprintf("%s & %d & %d & %.1f & %.1f & %.1f & %.1f & %.1f & %.1f \\\\", sh, g$r, g$n,
            g$rej_LRT, g$rej_ZC, g$rej_LL, g$rej_bonf, g$rej_simes, g$rej_fisher))
}
w("\\bottomrule\\end{tabular*}\\end{table}")

## ---- Table: power, the configurations discussed in the text ---------
LAB  <- c("2a"="2","2b"="2","3a"="3","3b"="3","4a"="4","4b"="4")
SETT <- c("2a"="(i)","2b"="(ii)","3a"="(i)","3b"="(ii)","4a"="(i)","4b"="(ii)")
pick <- data.frame(
  case = c("2a","2b","2b","3a","3b","4a","4b","4b"),
  r    = c(   5,   5,   4,   5,   5,   5,   5,   3),
  n    = c( 200, 200, 500, 200, 200, 200, 200,  50), stringsAsFactors = FALSE)
w("")
w("%%%% TABLE: power (excerpt)")
w("\\begin{table}[htbp]")
w("\\centering")
w("\\caption{Empirical power (\\%) at the nominal $5\\%$ level at the configurations",
  " discussed in the text; the full design is in Table~\\ref{Table_power_full} of",
  " Subsection~\\ref{app_power}}")
w("\\label{Table_power}")
w("\\begin{tabular*}{\\textwidth}{@{\\extracolsep\\fill}llrrrrrrr} \\toprule")
w("Case & setting & $r$ & $n$ & $G^{2}(H_{I})$ & $G^{2}(H_{ZC})$ & $G^{2}(H_{LL})$",
  " & Bonferroni & Fisher \\\\ \\midrule")
prevc <- ""
for (k in seq_len(nrow(pick))) {
  q <- pick[k, ]
  g <- main[main$case == q$case & main$shape == "skew" & main$r == q$r & main$n == q$n, ]
  c1 <- if (LAB[[q$case]] == prevc) "" else LAB[[q$case]]; prevc <- LAB[[q$case]]
  w(sprintf("%s & %s & %d & %d & %.1f & %.1f & %.1f & %.1f & %.1f \\\\",
            c1, SETT[[q$case]], q$r, q$n,
            g$rej_LRT, g$rej_ZC, g$rej_LL, g$rej_bonf, g$rej_fisher))
}
w("\\bottomrule\\end{tabular*}\\end{table}")

## ---- Table: power over the whole design (appendix, a longtable) -----
w("")
w("%%%% TABLE: power (full, appendix)")
w("\\begingroup")
w("\\setlength{\\LTleft}{0pt}\\setlength{\\LTright}{0pt}")
w("%% longtable fixes its caption width at 4in; widen it to the text block.")
w("\\setlength{\\LTcapwidth}{\\textwidth}")
w("\\begin{longtable}{@{\\extracolsep{\\fill}}llrrrrrrr}")
w("\\caption{Empirical power (\\%) at the nominal $5\\%$ level over the whole design,",
  " at the marginal labelled moderate in Table~\\ref{Table_size}}")
w("\\label{Table_power_full} \\\\")
HEAD <- paste("Case & setting & $r$ & $n$ & $G^{2}(H_{I})$ & $G^{2}(H_{ZC})$ &",
              "$G^{2}(H_{LL})$ & Bonferroni & Fisher \\\\")
w("\\toprule"); w(HEAD, " \\midrule"); w("\\endfirsthead")
w("\\multicolumn{9}{@{}l}{\\textit{Table~\\ref{Table_power_full} (continued)}} \\\\")
w("\\toprule"); w(HEAD, " \\midrule"); w("\\endhead")
w("\\midrule")
w("\\multicolumn{9}{r@{}}{\\textit{continued on the next page}} \\\\")
w("\\endfoot"); w("\\bottomrule"); w("\\endlastfoot")
prevc <- ""
for (cs in c("2a","2b","3a","3b","4a","4b")) {
  dd <- main[main$case == cs & main$shape == "skew", ]; dd <- dd[order(dd$r, dd$n), ]
  for (i in seq_len(nrow(dd))) {
    g <- dd[i, ]
    c1 <- if (LAB[[cs]] == prevc) "" else LAB[[cs]]; prevc <- LAB[[cs]]
    c2 <- if (i == 1) SETT[[cs]] else ""
    w(sprintf("%s & %s & %d & %d & %.1f & %.1f & %.1f & %.1f & %.1f \\\\",
              c1, c2, g$r, g$n, g$rej_LRT, g$rej_ZC, g$rej_LL, g$rej_bonf, g$rej_fisher))
  }
}
w("\\end{longtable}"); w("\\endgroup")

## ---- Table: sensitivity to the scores -------------------------------
midrank <- function(y, r) { n <- sum(y); p <- n_star(y, r) / (2 * n)
  cu <- cumsum(p); (c(0, cu[-r]) + cu) / 2 }
SETS <- function(y, r) {
  o <- list("$u_i = i$" = seq_len(r))
  if (r == 3) { o[["$(1,2,4)$"]] <- c(1,2,4); o[["$(1,3,4)$"]] <- c(1,3,4)
                o[["$(1,2,10)$"]] <- c(1,2,10)
  } else { o[["$(1,2,3,5)$"]] <- c(1,2,3,5); o[["$(1,3,4,5)$"]] <- c(1,3,4,5)
           o[["$(1,2,3,10)$"]] <- c(1,2,3,10) }
  o[["midranks"]] <- midrank(y, r); o[["normal"]] <- stats::qnorm(midrank(y, r)); o
}
DATA <- list("Table~\\ref{Table_3}" = list(y = c(4,17,8,21,20,6), r = 3),
             "Data 1" = list(y = c(14,85,194,65,81,272,136,263,318,42), r = 4),
             "Data 2" = list(y = c(45,85,74,33,65,93,69,62,52,24), r = 4))
w("")
w("%%%% TABLE: scores")
w("\\begin{table}[htbp]")
w("\\centering")
w("\\caption{The decomposition of the three data sets under six monotone scoring",
  " schemes, with $p$-values in parentheses}")
w("\\label{Table_scores}")
w("\\begin{tabular*}{\\textwidth}{@{\\extracolsep\\fill}llrrrr} \\toprule")
w("data & scores & $G^{2}(H_{ZC})$ & $G^{2}(H_{LL})$ & $G^{2}(H_{I})$ & $\\hat{\\theta}$ \\\\ \\midrule")
for (nm in names(DATA)) {
  d <- DATA[[nm]]; first <- TRUE
  for (sc in names(SETS(d$y, d$r))) {
    o <- decompose(d$y, d$r, SETS(d$y, d$r)[[sc]])
    w(sprintf("%s & %s & %.3f\\,(%.3f) & %.3f\\,(%.3f) & %.3f\\,(%.3f) & %.3f \\\\",
              if (first) nm else "", sc, o$G2[["ZC"]], o$p[["ZC"]],
              o$G2[["LL"]], o$p[["LL"]], o$G2[["I"]], o$p[["I"]], o$theta))
    first <- FALSE
  }
  if (nm != "Data 2") w("\\midrule")
}
w("\\bottomrule\\end{tabular*}\\end{table}")
close(con)
cat(sprintf("five tables written to %s (%d lines)\n", OUT, length(readLines(OUT))))
