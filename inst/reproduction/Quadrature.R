#!/usr/bin/env Rscript
## Quadrature.R
## ---------------------------------------------------------------------------
## The error of each evaluator of the logistic-normal integral against the
## variance argument v, and the time each takes.
##
## For each v, the largest absolute error over 801 values of eta in [-10, 10]
## of fixed-node Gauss-Hermite rules with 10 to 100 nodes, of the certified
## mixtures at K = 8, 11, and 22, and of two one-term approximations. The
## reference is adaptive Gauss-Kronrod quadrature, and the first column gives
## its largest difference from a 500-node Gauss-Hermite rule.
##
## The rules compete on equal terms, each being a fixed weighted sum of a
## sigmoid or of a normal distribution function, and the 11-node rule costs
## what the K = 11 mixture costs in function evaluations. Gauss-Hermite quadrature
## converges faster than any power of the number of nodes at fixed v, and its
## error grows with v, since the poles of the integrand approach the real axis
## in the quadrature variable. The mixture holds one bound at every v.
##
## Usage: Rscript Quadrature.R
## Requires statmod.
##
## R style: 4-space indentation, Map/Reduce/Filter instead of for-loops,
## `[[` instead of `$`, base pipe `|>`, functional replacement forms.

Sys.setenv(
    OPENBLAS_NUM_THREADS = "1",
    OMP_NUM_THREADS      = "1",
    MKL_NUM_THREADS      = "1"
)

find_pkg_root <- function() {
    candidates <- c(".", "..", file.path("..", ".."))
    hits <- Filter(function(p) file.exists(file.path(p, "DESCRIPTION")), candidates)
    if (length(hits) == 0L) NA_character_ else normalizePath(hits[[1L]])
}
if (!requireNamespace("LogisNormIntegral", quietly = TRUE)) {
    if (!requireNamespace("devtools", quietly = TRUE)) {
        stop("Install LogisNormIntegral first or run from a clone with devtools.")
    }
    pkg_root <- find_pkg_root()
    if (is.na(pkg_root)) {
        stop("Cannot find DESCRIPTION. Install LogisNormIntegral or run from the repo.")
    }
    devtools::load_all(pkg_root)
}
if (!requireNamespace("statmod", quietly = TRUE)) stop("Install the package statmod.")

## Fixed-node Gauss-Hermite quadrature against the standard normal law.
gh_eval <- function(n) {
    rule <- statmod::gauss.quad.prob(n, dist = "normal")
    function(eta, v) {
        Reduce(function(acc, j) {
            acc + rule[["weights"]][[j]] * plogis(eta + sqrt(v) * rule[["nodes"]][[j]])
        }, seq_along(rule[["nodes"]]), 0)
    }
}
reference <- function(eta, v) {
    mapply(function(e) {
        integrate(function(m) plogis(e + sqrt(v) * m) * dnorm(m), -Inf, Inf, rel.tol = 1e-13,
                  abs.tol = 1e-15, subdivisions = 500L)[["value"]]
    }, eta)
}
mixture <- function(K) {
    function(eta, v) LogisNormIntegral::logistic_normal_integral(eta, rep(v, length(eta)), K = K)
}
## Two one-term approximations. The probit form replaces the logistic function by
## pnorm(x / 1.6), as Gaynor, Schwartz, and Lin (2019) do in mediation. The attenuated
## logit keeps the logistic function and shrinks its argument, with
## c = 16 sqrt(3) / (15 pi), the formula of Zeger, Liang, and Albert (1988).
probit_one_term  <- function(eta, v) pnorm(eta / sqrt(1.6^2 + v))
attenuated_logit <- function(eta, v) plogis(eta / sqrt(1 + (16 * sqrt(3) / (15 * pi))^2 * v))
evaluators <- c(
    list(probit = probit_one_term, attenuated_logit = attenuated_logit),
    Map(gh_eval, c(GH_10 = 10L, GH_11 = 11L, GH_20 = 20L, GH_30 = 30L, GH_50 = 50L, GH_100 = 100L)),
    Map(mixture, c(mixture_K8 = 8L, mixture_K11 = 11L, mixture_K22 = 22L))
)
gh_500 <- gh_eval(500L)

eta_grid <- seq(-10, 10, length.out = 801L)
v_values <- c(0.01, 0.1, 0.5, 1, 2, 4, 9, 16, 25)
error_row <- function(v) {
    truth <- reference(eta_grid, v)
    c(v = v, reference = max(abs(truth - gh_500(eta_grid, v))),
      mapply(function(f) max(abs(f(eta_grid, v) - truth)), evaluators))
}
errors <- Map(error_row, v_values) |> Reduce(f = rbind)

cat("The error of each evaluator against the variance argument v\n")
cat("=========================================================================\n")
cat("  Largest absolute error over 801 values of eta in [-10, 10], against\n")
cat("  adaptive quadrature. 'reference' is the largest difference between that\n")
cat("  quadrature and a 500-node Gauss-Hermite rule.\n\n")
options(width = 200)
print(`rownames<-`(signif(errors, 2L), NULL))

## The smallest Gauss-Hermite rule whose error is at most the certified K = 11
## bound, for each v, searched over 5 to 400 nodes in steps of 5.
bound_11 <- LogisNormIntegral::minimax_error(11L)
nodes_needed <- function(v) {
    truth  <- reference(eta_grid, v)
    passes <- function(n) max(abs(gh_eval(n)(eta_grid, v) - truth)) <= bound_11
    hits   <- Filter(passes, seq(5L, 400L, by = 5L))
    if (length(hits) == 0L) NA_integer_ else hits[[1L]]
}
cat(sprintf("\nSmallest Gauss-Hermite rule, in steps of 5 nodes, with error at most %.3g:\n\n",
            bound_11))
print(data.frame(v = v_values, nodes = mapply(nodes_needed, v_values)), row.names = FALSE)

## Time for one million evaluations at v = 1, the median of five runs.
set.seed(6489L)
eta_time <- runif(1e6, -6, 6)
timing <- mapply(function(f) {
    median(mapply(function(i) system.time(f(eta_time, 1))[["elapsed"]], seq_len(5L)))
}, evaluators)
cat("\nSeconds for one million evaluations (median of five):\n\n")
print(round(timing, 2L))

invisible(errors)
