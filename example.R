# =====================================================================
#  Give it a table, get the result.
#
#  Run with:  Rscript example.R
# =====================================================================

source("R/intraclass.R")

## ---------------------------------------------------------------
## 1.  The simplest use: the upper triangle of the table, nothing else.
##     The number of categories is worked out from the length.
## ---------------------------------------------------------------
cat("\n##  Ishii (1960): progress of paired learning  ##\n\n")
ishii <- c(4, 17, 8,      # (1,1) (1,2) (1,3)
               21, 20,    #       (2,2) (2,3)
                    6)    #             (3,3)
print(analyse_intraclass(ishii))

## ---------------------------------------------------------------
## 2.  The same table as a matrix. The lower triangle may be zero or NA.
## ---------------------------------------------------------------
M <- matrix(NA, 3, 3)
M[1, ] <- c(4, 17, 8)
M[2, 2:3] <- c(21, 20)
M[3, 3] <- 6
stopifnot(all.equal(analyse_intraclass(M)$G2, analyse_intraclass(ishii)$G2))

## ---------------------------------------------------------------
## 3.  A full table of ordered pairs is folded automatically.
## ---------------------------------------------------------------
full <- matrix(c(4, 10, 5,
                 7, 21, 12,
                 3,  8,  6), 3, 3, byrow = TRUE)
stopifnot(all.equal(analyse_intraclass(full)$G2, analyse_intraclass(ishii)$G2))

## ---------------------------------------------------------------
## 4.  Other scores. G^2(H_I) does not involve them and does not move;
##     the two components into which it splits do.
## ---------------------------------------------------------------
cat("\n##  Sensitivity of the split to the scores  ##\n\n")
cat(sprintf("%-16s %9s %9s %9s\n", "scores", "G2(H_I)", "G2(H_ZC)", "G2(H_LL)"))
for (u in list(c(1, 2, 3), c(1, 2, 4), c(1, 3, 4), c(1, 2, 10))) {
  o <- analyse_intraclass(ishii, scores = u)
  cat(sprintf("%-16s %9.3f %9.3f %9.3f\n",
              paste0("(", paste(u, collapse = ","), ")"),
              o$G2["I"], o$G2["ZC"], o$G2["LL"]))
}

## ---------------------------------------------------------------
## 5.  The two artificial tables of the paper, showing the diagnostic use.
## ---------------------------------------------------------------
cat("\n##  Data 1: the LL association structure is at fault  ##\n\n")
print(analyse_intraclass(c(14, 85, 194, 65, 81, 272, 136, 263, 318, 42)))
cat("\n##  Data 2: the zero covariance structure is at fault  ##\n\n")
print(analyse_intraclass(c(45, 85, 74, 33, 65, 93, 69, 62, 52, 24)))
