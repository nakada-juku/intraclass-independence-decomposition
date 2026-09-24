# =====================================================================
#  Figure: the rate at which the discrepancy vanishes.
#
#  Companion to Fig_partition. That figure plots D itself, on which no
#  curve is straight; this one plots the standard deviation of D against
#  the sample size on both logarithmic scales, where the theoretical
#  rate sd(D) = O(n^{-1/2}) is a straight line of slope -1/2.
#
#  The point of the display is the departure from straightness: at the
#  marginals that have entered the asymptotic regime the points lie on a
#  line, and at the most concentrated one they do not.
#
#  Usage: Rscript make_figure_rate.R [outfile]
# =====================================================================

args <- commandArgs(trailingOnly = TRUE)
OUT  <- if (length(args) >= 1) args[1] else "../figures/Fig_rate.pdf"

## The three designs differ only in their sample sizes, over the same family of
## probability tables, so they are treated as one.
o <- do.call(rbind, lapply(
  Filter(file.exists, sprintf("results/%s.rds", c("order", "order2", "order3"))),
  readRDS))
## Sort by n, so that the lines are not drawn in the order the designs were bound
o <- o[order(o$shape, o$r, o$n), ]

LEV <- c("H1.00", "H0.95", "H0.85", "H0.75", "H0.50")
LAB <- sub("^H", "", LEV)
COL <- colorRampPalette(c("#1b7f79", "#2d4a7a", "#a5432a"))(length(LEV))
PCH <- c(16, 17, 15, 18, 4)
o  <- o[o$shape %in% LEV, ]
RS <- sort(unique(o$r)); NS <- sort(unique(o$n))

PT <- 9
if (grepl("\\.png$", OUT)) png(OUT, width = 6.5, height = 2.8, units = "in", res = 220, pointsize = PT) else
  pdf(OUT, width = 6.5, height = 2.8, pointsize = PT)
layout(matrix(1:3, 1, 3))
par(mar = c(3.4, 0.5, 1.8, 0.5), oma = c(0, 4.0, 0, 0.3), mgp = c(2.0, 0.55, 0), las = 1)

yt <- c(0.02, 0.05, 0.1, 0.2, 0.5, 1)
for (rr in RS) {
  d0 <- o[o$r == rr, ]
  plot(NA, xlim = range(NS) * c(0.85, 1.18), ylim = range(o$D_sd) * c(0.9, 1.1),
       log = "xy", axes = FALSE, xlab = "", ylab = "")
  abline(h = yt, col = "grey93")
  ## Reference line of slope -1/2, the theoretical order sd(D) = O(n^{-1/2})
  aref <- exp(mean(log(d0$D_sd[d0$n == 200])) + 0.5 * log(200))
  lines(NS, aref * NS^(-0.5), col = "grey45", lty = 2, lwd = 1.1)
  for (k in seq_along(LEV)) {
    d <- d0[d0$shape == LEV[k], ]
    if (!nrow(d)) next
    lines(d$n, d$D_sd, col = COL[k], lwd = 1.4)
    points(d$n, d$D_sd, col = COL[k], pch = PCH[k], cex = 0.6)
  }
  box(col = "grey55")
  at <- c(40, 100, 1000, 10000)
  axis(1, at = at, labels = c("40","100","1000",expression(10^4)),
       cex.axis = 0.76, tcl = -0.22)
  if (rr == RS[1]) {
    axis(2, at = yt, cex.axis = 0.82, tcl = -0.22)
    mtext(expression(paste("standard deviation of ", D)), side = 2, outer = TRUE,
          line = 2.7, las = 0, cex = 0.92)
  }
  mtext(bquote(r == .(rr)), side = 3, line = 0.4, cex = 0.95)
  mtext("sample size  n   (both scales logarithmic)", side = 1, line = 2.0, cex = 0.70)
  if (rr == RS[1])
    legend("bottomleft", bty = "n", cex = 0.72, title = expression(H(bold(p)) / log ~ r),
           legend = LAB, col = COL, lwd = 1.4, pch = PCH, pt.cex = 0.6,
           seg.len = 1.4, y.intersp = 0.95)
  if (rr == RS[length(RS)])
    legend("topright", bty = "n", cex = 0.72, legend = expression(slope ~ -1/2),
           col = "grey45", lty = 2, lwd = 1.1, seg.len = 1.8)
}
invisible(dev.off())
cat("figure written to", OUT, "\n")
