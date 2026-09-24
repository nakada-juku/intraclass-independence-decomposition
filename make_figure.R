# =====================================================================
#  Figure: the discrepancy of Theorem 2 against the sample size.
#
#  The quantity plotted is
#      D = G^2(H_I) - {G^2(H_ZC) + G^2(H_LL)}
#  itself, on a linear vertical scale. The horizontal scale is logarithmic
#  only because the design spans n = 40 to 300000; on a linear axis the
#  whole of the range in which D changes would occupy the left-hand
#  hundredth of each panel. It is not a linearising transformation: the
#  mean of D is O(n^{-1}) and its dispersion O(n^{-1/2}), so no curve here
#  is straight on any scale. The companion figure, which plots the
#  dispersion alone, is the one on which a power law does appear straight.
#
#  One panel per marginal and per r, so that the three need not be offset
#  against one another within a panel.
#
#  Usage: Rscript make_figure.R [outfile] [ymax]
# =====================================================================

args <- commandArgs(trailingOnly = TRUE)
OUT  <- if (length(args) >= 1) args[1] else "../figures/Fig_partition.pdf"
CAP  <- if (length(args) >= 2) as.numeric(args[2]) else 1.5

## The three designs differ only in their sample sizes, over the same family of
## probability tables, so they are treated as one.
o <- do.call(rbind, lapply(
  Filter(file.exists, sprintf("results/%s.rds", c("order", "order2", "order3"))),
  readRDS))
## Sort by n, so that the lines are not drawn in the order the designs were bound
o <- o[order(o$shape, o$r, o$n), ]

LEV <- c("H1.00", "H0.75", "H0.50")
LAB <- sub("^H", "", LEV)
COL <- c("#1b7f79", "#2d4a7a", "#a5432a")
o  <- o[o$shape %in% LEV, ]
RS <- sort(unique(o$r)); NS <- sort(unique(o$n))

TICK <- 0.5
lo <- floor(min(o$D_p05) / TICK) * TICK
hi <- CAP
tk <- seq(lo, hi, by = TICK)

PT <- 9
if (grepl("\\.png$", OUT))
  png(OUT, width = 6.5, height = 5.4, units = "in", res = 220, pointsize = PT) else
  pdf(OUT, width = 6.5, height = 5.4, pointsize = PT)
par(mfrow = c(3, 3), mar = c(0.4, 0.5, 0.4, 0.5),
    oma = c(3.8, 5.8, 2.2, 0.6), mgp = c(2.0, 0.55, 0), las = 1)

for (k in seq_along(LEV)) for (rr in RS) {
  d <- o[o$shape == LEV[k] & o$r == rr, ]
  plot(NA, xlim = range(NS) * c(0.75, 1.3), ylim = range(tk), log = "x",
       axes = FALSE, xlab = "", ylab = "")
  abline(h = tk, col = "grey93"); abline(h = 0, col = "grey35", lwd = 1.0)
  if (nrow(d)) {
    sp   <- diff(range(tk))
    over <- d$D_p95 > hi
    up   <- ifelse(over, hi - 0.055 * sp, d$D_p95)
    nv   <- !over
    if (any(nv))
      arrows(d$n[nv], d$D_p05[nv], d$n[nv], d$D_p95[nv], code = 3, angle = 90,
             length = 0.018, col = COL[k], lwd = 1.0)
    if (any(over)) {
      arrows(d$n[over], pmax(d$D_p05[over], lo), d$n[over], up[over],
             code = 1, angle = 90, length = 0.018, col = COL[k], lwd = 1.0)
      arrows(d$n[over], up[over], d$n[over], hi - 0.012 * sp,
             length = 0.05, angle = 20, col = COL[k], lwd = 1.0)
    }
    lines(d$n, d$D_mean, col = COL[k], lwd = 1.6)
    points(d$n, d$D_mean, col = COL[k], pch = 16, cex = 0.5)
  }
  box(col = "grey55")
  if (k == length(LEV)) {
    ## n spans more than two decades, so the ticks are placed at powers of ten
    at <- c(40, 100, 1000, 10000)
    axis(1, at = at, labels = c("40", "100", "1000", expression(10^4)),
         cex.axis = 0.76, tcl = -0.22)
  }
  if (rr == RS[1]) {
    axis(2, at = tk, cex.axis = 0.78, tcl = -0.22)

    mtext(bquote(H(bold(p)) / log ~ r %~~% .(LAB[k])), side = 2, line = 2.5,
          las = 0, cex = 0.72)
  }
  if (k == 1) mtext(bquote(r == .(rr)), side = 3, line = 0.5, cex = 0.82)
}
mtext("sample size  n   (logarithmic)", side = 1, outer = TRUE, line = 2.2, cex = 0.78)
mtext(expression(D == G^2 * (H[I]) - group("{", G^2 * (H[ZC]) + G^2 * (H[LL]), "}")),
      side = 2, outer = TRUE, line = 4.2, las = 0, cex = 0.86)
invisible(dev.off())
cat("figure written to", OUT, "\n")
