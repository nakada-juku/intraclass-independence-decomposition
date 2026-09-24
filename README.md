# Decomposition of the independence model for two-way intraclass contingency tables

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22928660.svg)](https://doi.org/10.5281/zenodo.22928660)

R code accompanying *Decomposition of Independence Model for Two-Way Intraclass
Contingency Tables with Ordered Categories* (Nakada, Tahata, Tomizawa and Shinoda).

It fits the independence model `H_I` of Ishii (1960), the zero covariance model
`H_ZC` and the linear-by-linear association model `H_LL` to an intraclass table,
and reports the likelihood ratio tests together with the partition of Theorem 2.

Only base R and `stats` are required. `Rsolnp` is *not* needed: the constrained
fits are solved directly (see "How the models are fitted" below).

## Quick start

```r
source("R/intraclass.R")

## Ishii (1960): progress of paired learning, r = 3, n = 76 pairs
analyse_intraclass(c(4, 17, 8, 21, 20, 6))
```

```
Intraclass contingency table:  r = 3,  n = 76
Scores: 1, 2, 3

                           G2 df p.value
H_I  (independence)     0.243  3   0.970
H_ZC (zero covariance)  0.190  1   0.663
H_LL (linear-by-linear) 0.049  2   0.976

theta.hat = 1.11
G2(H_I) - {G2(H_ZC) + G2(H_LL)} = 0.0046

Combined p-values for H_I from the two components:
  Bonferroni 1.000   Simes 0.976   Fisher 0.929
```

`Rscript example.R` runs this and several variations.

## Giving it the data

`analyse_intraclass(y, scores = NULL)` accepts the table in any of three forms
and **works out the number of categories itself**:

| form | example |
|---|---|
| upper triangle, including the diagonal, in row order | `c(4, 17, 8, 21, 20, 6)` |
| an `r × r` matrix with the counts in its upper triangle (lower triangle `0` or `NA`) | see `example.R` |
| a full `r × r` table of ordered pairs, which is folded automatically | see `example.R` |

`scores` takes any strictly increasing vector of length `r`; it defaults to
`1, ..., r`. The statistic `G2(H_I)` does not involve the scores and does not
change with them, but the two components into which Theorem 2 splits it do, so
it is worth trying more than one scoring.

## Files

```
R/intraclass.R    the three models, their tests, and analyse_intraclass()
R/populations.R   the probability tables used in the simulation study
R/design.R        the designs of the study and the code for one condition
run_simulation.R  runner: parallel, resumable, writes each condition as it finishes
collect_results.R assembles the per-condition files into one table
make_tables.R     every LaTeX table of the paper
make_figure.R     Figure 2, the discrepancy against the sample size
make_figure_rate.R Figure 3, its dispersion on both logarithmic scales
verify_paper.R    checks the manuscript against the results
scores.R          the sensitivity analysis of Section 5
example.R         the examples above
tests.R           the checks described below
```

## Reproducing the simulation study

Four designs are used. `main` is the study of size and power; `order`,
`order2` and `order3` together make up the study of the partition, the second
and third adding sample sizes to the first.

```sh
Rscript run_simulation.R main   5000 4     # design, replications, workers
Rscript run_simulation.R order  5000 4
Rscript run_simulation.R order2 5000 4
Rscript run_simulation.R order3 5000 4
for d in main order order2 order3; do Rscript collect_results.R $d; done
Rscript make_tables.R && Rscript make_figure.R && Rscript make_figure_rate.R
Rscript verify_paper.R
```

Each condition is written to `results/<design>/cond_XXXX.rds` **as soon as it
finishes**, and a condition whose file already exists is skipped. The run can
therefore be interrupted and restarted at any point without losing work, and
nothing is accumulated in memory across conditions. On a laptop it is worth
starting it under `nice -n 10`. The default number of workers is half the
available cores, so the machine stays usable. A subset may be given as a
fourth argument, as in `Rscript run_simulation.R order 5000 4 85-105`.

Every replication is seeded individually, from the design and the condition
number, so results do not depend on the number of workers. **A consequence is
that the sample sizes of a design must not be altered in place**: inserting one
shifts the condition numbers, and the stored results would then belong to a
different random stream than the code would generate. That is why the extra
sample sizes are separate designs rather than additions to `order`, and why
`collect_results.R` refuses to assemble files that do not match the grid.

The last two designs take the longest, not because of the sample size but
because of the empty cells: the barrier iterations in the zero covariance fit
are what cost, and they disappear as the expected frequencies grow.

## How the models are fitted

`H_I` has a closed form, `alpha_i = (n_{i+} + n_{+i}) / (2n)`.

`H_LL` is fitted by unconstrained maximisation over `(log alpha, log theta)`,
the normalisation being imposed by construction.

`H_ZC` needs more care. The constraint `Cov(U,V) = g'pi - (s'pi)^2 = 0` is
quadratic and not convex, so a penalty on it can jump between branches and
return a point that is not the maximum. It becomes linear once the mean
`m = s'pi` is held fixed, because then

```
Cov = 0   <=>   s'pi = m   and   g'pi = m^2 ,
```

so the fit is obtained by maximising a concave log-likelihood over the simplex
subject to two linear constraints -- solved in its three-dimensional convex
dual -- and then maximising the profile over the scalar `m`. A cell with no
observations may still carry probability under this model, because it relaxes
the two moment constraints, so the off-support cells are handled by a
logarithmic barrier rather than being fixed at zero.

## Checks

`Rscript tests.R` verifies that the code reproduces the values printed in the
paper and that the fits satisfy the relations they must satisfy.
