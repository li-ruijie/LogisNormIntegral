# LogisNormIntegral

<!-- badges: start -->
[![R-CMD-check](https://github.com/li-ruijie/LogisNormIntegral/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/li-ruijie/LogisNormIntegral/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

`LogisNormIntegral` is an R package for evaluating the **logistic-normal integral**, the
expected logistic response of a Gaussian latent variable, with a proved bound on the error
of every evaluation.

## The problem

Averaging a logistic response over a normally distributed linear predictor gives

$$
I(\eta, v) \;=\; \mathbb{E}\big[\Lambda(Y)\big], \qquad Y \sim N(\eta, v), \qquad
\Lambda(y) = \frac{1}{1 + e^{-y}} .
$$

This integral has no exact closed form, so it is ordinarily evaluated by adaptive
quadrature or by inner Monte Carlo sampling. Both are costly in computation time when the integral sits
inside an outer loop over parameter draws, which is where applications place it. The
inexpensive alternative in common use, a second-order Taylor expansion in $v$, carries a
bias that does not vanish as the outer sample grows.

This package supplies a third route. A fixed K-term normal-scale mixture collapses the
integral to a finite sum of normal distribution functions, evaluable in closed form, whose
absolute error is bounded uniformly over all $(\eta, v)$ by one proved constant
Delta_K\*. That bound is a property of the shipped constants rather than an estimate
recomputed at each call.

## What the package provides

- **Certified minimax normal-scale mixtures.** A K-term mixture
  $\sum_i w_i \Phi(s_i z)$ for the logistic CDF, with uniform absolute error below a
  proved bound Delta_K\*. For K = 1, ..., 11 the shipped weights and scales are the
  certified global minimax, proved to equioscillate at 2K points and to be the unique
  minimiser via varisolvency. Relaxed Bernstein-Akhiezer certificates at K = 14 to 17
  reach Delta = 7.2e-12 to 2.5e-12 with all weights positive, and a signed
  K = 22 frontier reaches
  Delta = 3.6e-15, within a factor of about thirty-three of the IEEE double-precision
  unit roundoff 2^-53.
- **Closed-form asymptotic series.** Small- and large-variance expansions from the
  exact duality between the integral and Zwegers' Mordell integral.

## Application: causal mediation

The integral is the computational kernel of causal mediation analysis with a binary
outcome and a continuous mediator. The expected potential outcomes, and hence the natural
indirect effect, are logistic-normal integrals in the fitted model parameters, so a
quasi-Bayesian analysis evaluates one per outer draw. A certified mixture evaluation
replaces the inner quadrature or sampling, and its error is bounded by Delta_K\* rather
than by a sampling rate. Against inner Monte Carlo at $M_{\rm inner} = 10^4$,
`inst/reproduction/Example1.R` measures a floor of 380 times faster at the K = 22 default.
The second-order Taylor expansion of Valeri and VanderWeele is comparable in speed and
carries a bias at each of the three outcome prevalences that script tests, flipping sign
across them and largest at high prevalence. `predict_indirect_effect_qb()` is the entry
point, shown under Usage below.

## Installation

The package needs only base R (>= 3.5.0). Install a tagged source release straight from
GitHub:

```r
# install.packages("remotes")
remotes::install_github("li-ruijie/LogisNormIntegral@v0.2.5")
```

Every [release](https://github.com/li-ruijie/LogisNormIntegral/releases) also carries a
pre-built source tarball, installable without any build tools:

```r
install.packages(
  "https://github.com/li-ruijie/LogisNormIntegral/releases/download/v0.2.5/LogisNormIntegral_0.2.5.tar.gz",
  repos = NULL, type = "source")
```

## Usage

```r
library(LogisNormIntegral)

# The logistic-normal integral E[Lambda(Y)], Y ~ N(eta, v)
logistic_normal_integral(eta = 0, v = 1)             # 0.5 by symmetry
logistic_normal_integral(c(-1, 0, 1), v = 2, K = 11)

# Certified error bound, weights, and scales for a given order
minimax_error(11)                                    # 1.26e-11
minimax_weights(8)
minimax_scales(8)

# The logistic-CDF approximation itself
phi_minimax_mixture(2, K = 10)                       # close to plogis(2)

# Closed-form asymptotic series, branch selected per argument
phi_analytic(0.5, v = 1,  order = 2)                 # small-v series
phi_analytic(10,  v = 20, order = 2)                 # large-v continuation
```

By default the mixture functions use the signed K = 22 frontier (the lowest certified
error). Pass `K` to trade accuracy for fewer atoms. The accepted orders are
`1:11` (strict), `14` to `17` (relaxed BA, all weights positive), and `22`
(signed). Choose one of 14 to 17 when a strictly positive mixture is wanted
below the K = 11 error.

### Fast quasi-Bayesian causal mediation

`predict_indirect_effect_qb()` returns the natural indirect effect (NIE) and expected
potential outcomes for a binary outcome and continuous mediator, with properly
asymmetric quasi-Bayesian intervals:

```r
set.seed(42)
n     <- 500
treat <- rbinom(n, 1, 0.5)
covar <- rnorm(n)
med   <- 0.5 + 1.2 * treat + 0.3 * covar + rnorm(n, sd = 0.8)
y     <- rbinom(n, 1, plogis(-1 + 0.5 * treat + 1.5 * med - 0.2 * covar))
df    <- data.frame(y, treat, med, covar)

mod_m <- lm(med ~ treat + covar, data = df)
mod_y <- glm(y ~ treat + med + covar, data = df, family = binomial)

res <- predict_indirect_effect_qb(
  outcome_model  = mod_y,
  mediator_model = mod_m,
  covar_vals     = data.frame(covar = 0),
  treatment_var  = "treat",
  mediator_var   = "med")
res$NIE
```

## Certified accuracy

Certified uniform absolute error over the whole real line. For K = 1, ..., 11 this is the
strict global minimax Delta_K\*, proved to equioscillate at 2K points and to be the unique
minimiser. Relaxed Bernstein-Akhiezer certificates at K = 14 to 17, omitted from the table,
carry audited bounds from 7.3e-12 down to 2.6e-12 with all weights positive. The K = 22
entry is the signed relaxed frontier, a certified bound on that particular mixture rather
than a global minimax:

| K            | 1      | 2      | 4      | 6      | 8      | 10      | 11      | 22 (signed) |
|--------------|--------|--------|--------|--------|--------|---------|---------|-------------|
| Certified Δ  | 9.5e-3 | 5.1e-4 | 4.7e-6 | 8.4e-8 | 2.1e-9 | 6.6e-11 | 1.3e-11 | 3.6e-15     |

Displayed values are rounded to two significant figures, so five of the eight round down and are
not themselves valid upper bounds. Use `minimax_error(K)` for the value to double precision, and
`inst/certs/deltas.csv` for the audited supremum at its certified precision, which is the figure
to quote when a guarantee is required.

The mixtures ship twice, as the LazyData object under `data/` that `minimax_weights()` and
`minimax_scales()` read, and as a plain-text audit deposit under `inst/certs/`. `minimax_error()`
returns the certified bound itself. `inst/certs/README.md` documents the provenance of both copies
and carries a verification snippet.

## Reproduction

`inst/reproduction/` regenerates the certified tables and figures from the shipped cert data:

- `Table1.R`: the certified minimax constants for K = 1, ..., 11.
- `Frontier.R`: the certified frontier over all three families, strict K = 1, ..., 11,
  relaxed BA at K = 14 to 17, and the signed K = 22 bound.
- `Example1.R`: quasi-Bayesian inner-integral speed, bias, and accuracy comparison.
- `Coverage.R`: repeated-sampling coverage of the quasi-Bayesian NIE intervals.
- `PirjolTiming.R`: per-point accuracy and cost of the corrected Pirjol (2013)
  Poisson evaluator against the certified K = 11 mixture (needs a clone of the
  upstream `dan-pirjol/logisticNormal` package, see `inst/reproduction/README.md`).

## Mathematical background

The package vignette derives the two families of approximations, states the certification
standard the shipped constants meet, and sets out the equioscillation and global-uniqueness
results:

```r
vignette("certified-minimax", package = "LogisNormIntegral")
```

Its source is `vignettes/certified-minimax.Rmd`.

## License

GNU Affero General Public License v3.0 (`LICENSE.md` in the source repository).
Copyright (c) 2026 Ruijie Li.
