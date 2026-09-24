# =====================================================================
#  Assemble the per-condition files written by run_simulation.R into a
#  single table, and report anything that did not finish.
#
#  Usage: Rscript collect_results.R [design]
# =====================================================================

source("R/intraclass.R")
source("R/populations.R")
source("R/design.R")

args   <- commandArgs(trailingOnly = TRUE)
DESIGN <- if (length(args) >= 1) args[1] else "main"
OUTDIR <- file.path("results", DESIGN)
GRID   <- make_grid(DESIGN)

files <- sprintf("%s/cond_%04d.rds", OUTDIR, seq_len(nrow(GRID)))
have  <- file.exists(files)
if (!any(have)) stop("no results in ", OUTDIR)

res <- do.call(rbind, lapply(files[have], readRDS))

## Check that each stored file belongs to the grid row it is named for. The
## condition number fixes the random stream through set.seed(10000 * i + b), so
## a file that has merely been renumbered holds results its seed would not
## reproduce.
chk <- GRID[have, c("case", "r", "n", "shape")]
got <- res[, c("case", "r", "n", "shape")]
mism <- which(!Reduce(`&`, lapply(names(chk), function(k)
  as.character(chk[[k]]) == as.character(got[[k]]))))
if (length(mism))
  stop(sprintf("the files for conditions %s do not match the grid; they look renumbered",
               paste(which(have)[mism], collapse = ", ")))
out <- file.path("results", paste0(DESIGN, ".rds"))
saveRDS(res, out)

cat(sprintf("%d/%d conditions collected into %s\n", sum(have), nrow(GRID), out))
if (any(!have))
  cat("conditions still to do:", paste(which(!have), collapse = ", "), "\n")
bad <- res$fail_rate > 0 | res$invalid_rate > 0
if (any(bad)) {
  cat("\nconditions in which a fit gave trouble:\n")
  print(res[bad, c("case", "r", "n", "shape", "fail_rate", "invalid_rate", "nonconv_rate")])
} else cat("no failed or invalid fits\n")
