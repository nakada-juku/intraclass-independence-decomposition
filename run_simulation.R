# =====================================================================
#  Runner for the simulation study.
#
#  Design notes:
#    * conditions are run in parallel, one worker per condition;
#    * each condition writes its own small file as soon as it finishes,
#      and the worker returns NULL, so nothing accumulates in the parent;
#    * a condition whose file already exists is skipped, so the run is
#      resumable and can be interrupted at any time without loss;
#    * `collect_results.R` assembles the per-condition files afterwards.
#
#  Usage
#    Rscript run_simulation.R [design] [B] [cores] [subset]
#      design : "main" (default), "curve" or "order"
#      B      : replications per condition (default 5000)
#      cores  : workers (default: half the cores, to leave the machine usable)
#      subset : optional condition numbers to restrict the run to, given as
#               a comma separated list of numbers and ranges, e.g. "106-126"
#               or "1,4,7-9". Conditions outside the subset are left for a
#               later call. Useful for running the conditions one needs first
#               when the whole design would take longer than one is willing
#               to wait.
#
#  On a laptop it is worth starting this under `nice`:
#      nice -n 10 Rscript run_simulation.R main 5000 4
# =====================================================================

suppressMessages(library(parallel))
source("R/intraclass.R")
source("R/populations.R")
source("R/design.R")

args   <- commandArgs(trailingOnly = TRUE)
DESIGN <- if (length(args) >= 1) args[1] else "main"
B      <- if (length(args) >= 2) as.integer(args[2]) else 5000
CORES  <- if (length(args) >= 3) as.integer(args[3]) else max(1, floor(detectCores() / 2))

## A separate stream per design. The existing designs keep 0, so their results
## are unchanged.
SEEDBASE <- switch(DESIGN, order2 = 5000000L, order3 = 9000000L, 0L)
GRID   <- make_grid(DESIGN)
OUTDIR <- file.path("results", DESIGN)
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

fname <- function(i) file.path(OUTDIR, sprintf("cond_%04d.rds", i))
todo  <- which(!file.exists(vapply(seq_len(nrow(GRID)), fname, character(1))))

## Expand "106-126" or "1,4,7-9" into a set of condition numbers
parse_subset <- function(s) {
  parts <- strsplit(s, ",", fixed = TRUE)[[1]]
  unlist(lapply(parts, function(p) {
    p <- trimws(p)
    if (grepl("^[0-9]+-[0-9]+$", p)) {
      ab <- as.integer(strsplit(p, "-", fixed = TRUE)[[1]]); seq(ab[1], ab[2])
    } else if (grepl("^[0-9]+$", p)) {
      as.integer(p)
    } else stop("malformed subset: ", p)
  }))
}
SUBSET <- if (length(args) >= 4) parse_subset(args[4]) else NULL
if (!is.null(SUBSET)) {
  bad <- setdiff(SUBSET, seq_len(nrow(GRID)))
  if (length(bad)) stop("condition numbers out of range: ", paste(bad, collapse = ", "))
  todo <- intersect(todo, SUBSET)
}

cat(sprintf("design = %s / %d conditions (%d to do) / B = %d / workers = %d%s\n",
            DESIGN, nrow(GRID), length(todo), B, CORES,
            if (is.null(SUBSET)) "" else sprintf(" / subset = %s", args[4])))
if (!length(todo)) { cat("nothing to do\n"); quit(save = "no") }

t0 <- Sys.time()
invisible(mclapply(todo, function(i) {
  res <- tryCatch(run_condition(i, GRID, B, seed_base = SEEDBASE),
                  error = function(e) NULL)
  if (!is.null(res)) {
    saveRDS(res, fname(i))            # write it out the moment it is finished
    cat(sprintf("  [%s] condition %d/%d done (%.1f min)\n", format(Sys.time(), "%H:%M:%S"),
                i, nrow(GRID), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
    flush.console()
  }
  NULL                                 # return nothing, so the parent holds nothing
}, mc.cores = CORES, mc.preschedule = FALSE))

done <- sum(file.exists(vapply(seq_len(nrow(GRID)), fname, character(1))))
cat(sprintf("all done, %d/%d conditions (%.1f min)\n", done, nrow(GRID),
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
