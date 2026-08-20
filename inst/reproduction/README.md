# Manuscript reproduction scripts

Each script regenerates one artifact (table or figure) of the accompanying
manuscript from the cert data shipped with the package.

## Scripts

- `Table1.R`: the manuscript's certified-constants table. Certified minimax constants
  $\Delta_K^\star$ for $K = 1, \dots, 11$, with the closure pipeline
  for each $K$.
- `Frontier.R`: the certified frontier over all three shipped families, the
  strict global minimax for $K = 1, \dots, 11$, the relaxed Bernstein-Akhiezer
  certificates at $K = 14$ to $K = 17$, and the signed $K = 22$ bound
  ($\Delta = 3.604 \times 10^{-15}$), with an empirical max-error check at
  $K = 22$.
- `Example1.R`: the application-section worked example and its method-comparison table, a quasi-Bayesian
  speed/bias/accuracy comparison of inner-integral methods (minimax at the
  $K = 22$ default and at $K = 17$ and $K = 11$, V&VW Taylor, inner Monte
  Carlo, and optionally Pirjol's pointwise Poisson evaluator at truncation
  orders 10 and 30) across three outcome-prevalence regimes, against an
  adaptive-quadrature ground truth. Output: `Example1.output.txt`.
  The Pirjol rows need the upstream sources, which are third-party and not
  shipped. Clone them and pass the R directory to include those rows:
  ```sh
  git clone https://github.com/dan-pirjol/logisticNormal
  Rscript Example1.R logisticNormal/R
  ```
  Without a path the sweep runs without that engine and says so. The
  loop-bound fix of pull request 1 is applied in memory, so the upstream
  files are used unmodified on disk.
- `Coverage.R`: the application-section repeated-sampling claim, frequentist coverage
  of the nominal 95% quasi-Bayesian NIE intervals over 2000 simulated data sets
  per outcome-prevalence regime, for the certified minimax $K = 11$ mixture and,
  off the same outer draws, the V&VW Taylor expansion. Output:
  `Coverage.output.txt`.
- `PirjolTiming.R`: the application-section pointwise comparison. Per-point accuracy and
  per-call cost of the corrected Pirjol (2013) Poisson evaluator at truncation
  orders 30 and 10 against adaptive Gauss-Kronrod reference values, and of the
  certified minimax $K = 11$ sum, on the eight-point grid of the upstream bug
  report. Output: `PirjolTiming.output.txt`.

## Running

After `R CMD INSTALL ../..` or `devtools::install("../..")`:

```sh
Rscript Table1.R
Rscript Frontier.R
Rscript Example1.R
Rscript Coverage.R
```

`PirjolTiming.R` additionally needs a clone of the upstream
[`dan-pirjol/logisticNormal`](https://github.com/dan-pirjol/logisticNormal)
package, which is not shipped with this repository:

```sh
git clone https://github.com/dan-pirjol/logisticNormal
Rscript PirjolTiming.R logisticNormal/R
```

The script applies the one-line loop-bound fix of upstream pull request 1 in
memory, so the cloned files are used unmodified on disk.

Or directly from a clone without installing:

```sh
Rscript -e 'devtools::load_all("../.."); source("Table1.R")'
```

`Table1.R`, `Frontier.R`, `Example1.R`, and `Coverage.R` all support the same
two paths to the package, the installed `LogisNormIntegral` or
`devtools::load_all` from a clone, and stop with the same message when neither
is available. Devtools is therefore needed only for the clone path. `Table1.R`,
`Frontier.R`, and `Coverage.R` read the certified constants via package
accessors (`minimax_error(K)`, `minimax_weights(K)`, ...), while `Example1.R`
reaches them through `LogisNormIntegral::logistic_normal_integral()`. `Table1.R`
and `Frontier.R` print their tables to stdout in a near-publication format.
`PirjolTiming.R` loads no package at all, reading `../certs/minimax_k11.csv`
directly.

## Software environment

The package depends only on base R (>= 3.5.0). `Table1.R`, `Frontier.R`, and
`Coverage.R` use only base R and `LogisNormIntegral`. `Example1.R` additionally
needs `statmod` (Gauss-Hermite cross-check of the ground-truth integral), and
optionally the cloned upstream `logisticNormal` sources for its pointwise rows.
`PirjolTiming.R` uses base R, the shipped $K = 11$ certificate, and the cloned
upstream sources named above. Running the four package-using scripts from a
clone rather than against an installed package also needs `devtools`. The
multivariate-normal parameter draws are made with base `chol()`, so no `MASS`
dependency remains anywhere in the package.

## Expected output

`Table1.R` prints 11 rows. The first 8 reproduce Monahan & Stefanski 1992
Table 18.4.1 to 5 significant figures. The last 3 report the certifications of
this work to 16 significant figures, the ceiling a double carries, with the
audited full-precision bounds in `../certs/deltas.csv`. The monotonic ratio
check at the bottom shows ratios increasing from 0.054 (K = 1→2) to 0.189
(K = 10→11).

`Frontier.R` prints 16 rows (K = 1..11 strict, K = 14 to 17 relaxed BA, and
K = 22 signed) and the K = 22 empirical fp64 max-error (~4e-15, the certified
bound plus fp64 evaluation noise).

`Example1.R` prints a pointwise-accuracy table (K = 8/11/17 vs Monte Carlo) and a
quasi-Bayesian NIE comparison across three outcome-prevalence regimes. The minimax
$K = 11$ NIE is ~0.41 / 0.32 / 0.24 (low / moderate / high) and agrees with $K = 17$
to >= 10 decimals (the cert hierarchy is below QB sampling resolution), while V&VW
Taylor is biased by -13 / +9 / +25 pp. When the Pirjol engine is enabled it
reproduces the mixture to all six printed decimals in the low and moderate
regimes, and is undefined on 3 of 10000 draws at order 10 and 4 at order 30 in
the high regime. That evaluator returns a non-finite value at both ends of its
domain. The lower edge is at $v$ below about $0.0265$ at either order, is
order-independent, and sits where the total mediator coefficient is near zero.
All three order-10 draws and three of the four order-30 draws fall below it.
Both upper edges depend on the mean of the inner normal as well as on $v$. The
order-10 edge lies near $v = 142$ at a zero mean, above every draw of this
regime, so it is not exercised here. At order 30 the upper edge lies within
the sampled range, so no single threshold in $v$ describes it.
The fourth order-30 draw is the largest sampled variance of the regime,
$v = 56.5$, the only draw above the lower edge on which the order-30 evaluator
failed. Those draws are dropped and the count is printed. Needs `statmod`.

`Coverage.R` prints a 3-row coverage table. The certified minimax $K = 11$ QB
interval covers the fixed population NIE at close to the nominal 95% in all three
regimes (95.8 / 95.7 / 95.2%, low / moderate / high, as captured), while the V&VW
Taylor coverage is 4.9 / 99.2 / 92.2% off the same outer draws. Since the draws
are formed from a Cholesky factor the percentages reproduce across BLAS builds
(see the platform note in `Example1.R`). Runs in about 80 s.

`PirjolTiming.R` prints an eight-row per-point table and a per-call cost
summary. In the captured output the corrected Poisson evaluator at order 30
matches the Gauss-Kronrod reference to at most 2.3e-16 at about 21.5 us per
call, order 10 to at most 1.6e-12 at about 8.7 us, and the certified $K = 11$
mixture to at most 2.4e-14, within its certified bound, at about 1.7 us.
Timings vary with hardware.
