## Coverage.R: Frequentist (repeated-sampling) coverage of the quasi-Bayesian
##             NIE intervals built on the certified minimax mixture.
##
## Run from the package root or from inst/reproduction/:
##   Rscript inst/reproduction/Coverage.R
##
## Captured output: see Coverage.output.txt next to this script.
##
## Platform note. Seeded, and the parameter draws use `rmvnorm_chol` below,
## a Cholesky factorisation of the covariance, which is unique for a positive
## definite matrix and so reproduces across BLAS and LAPACK builds. This replaced
## MASS::mvrnorm on 2026-08-05, whose eigen() decomposition
## returned eigenvector signs that differ between builds. Coverage was the less
## exposed of the two scripts in either case, a rate being an average over 2000
## replicates. Measured against the earlier capture on a second platform, the
## largest move was 0.2 percentage points, inside the binomial Monte Carlo error
## the script prints alongside each rate. See the fuller note in Example1.R.
##
## ---------------------------------------------------------------------------
## What this script demonstrates
## ---------------------------------------------------------------------------
##
## Example1.R shows that the certified mixture reproduces the exact data-
## conditional NIE and yields asymmetric quasi-Bayesian (QB) intervals. Those
## intervals are credible intervals, valid in the Bayesian sense under the
## Bernstein--von Mises correspondence. This script asks the complementary
## frequentist question: across repeated samples from a known data-generating
## process, how often does the nominal 95% QB interval contain the fixed
## population NIE?
##
## For each of the three outcome-prevalence regimes of Example1.R
## (theta_0 = -3, -1, 0), whose marginal P(Y=1) is ~ 0.38, 0.69, 0.82 on the
## single data set of that script and 0.40, 0.67, 0.79 averaged over the
## replicates here, we
##   1. draw a fresh data set of n = 500 (treat -> mediator -> binary outcome),
##   2. fit the outcome (logistic, with treat:mediator interaction) and
##      mediator (linear) models by maximum likelihood,
##   3. build the 95% QB interval for the NIE (R_qb outer parameter draws,
##      inner integral by the certified mixture and, for contrast, by the
##      Valeri--VanderWeele second-order Taylor expansion),
##   4. record whether each interval covers the population NIE.
## Empirical coverage over N_sim replicates, with its binomial Monte Carlo
## standard error, is reported per regime and method.
##
## The inner integral for the mixture reproduces
## LogisNormIntegral::logistic_normal_integral(eta, v, K = 11) via the
## closed-form scale-mixture convolution
##   I(eta, v) = E_M[Lambda(eta + b M)]
##             = sum_i w_i * Phi(s_i * eta / sqrt(1 + s_i^2 * v)),   M ~ N(., v/b^2),
## using the certified strict K = 11 constants from the package accessors
## minimax_weights(11L) / minimax_scales(11L). Using K = 22 (the package default)
## changes no reported figure, the inner-integral error being bounded by
## Delta_K^* <= 1.3e-11, far below QB sampling resolution.
##
## R style: 4-space indentation, Map/Reduce/Filter instead of for-loops,
## `[[` instead of `$`, base pipe `|>`, pure functions, functional
## replacement forms.

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
## Sign-determinate multivariate normal draw. See the platform note above and the
## matching helper in Example1.R. The Cholesky factor of a positive definite matrix
## is unique, so the draws reproduce across BLAS and LAPACK builds, where the
## eigendecomposition MASS::mvrnorm uses does not. The pivoted branch preserves the
## tolerance of a numerically singular covariance and is also sign-determinate.
rmvnorm_chol <- function(n, mu, Sigma) {
    p <- length(mu)
    R <- tryCatch(chol(Sigma), error = function(e) NULL)
    if (is.null(R)) {
        pc  <- chol(Sigma, pivot = TRUE)
        R   <- pc[, order(attr(pc, "pivot")), drop = FALSE]
    }
    out <- matrix(stats::rnorm(n * p), nrow = n, ncol = p) %*% R
    out <- sweep(out, 2L, mu, "+")
    `dimnames<-`(out, list(NULL, names(mu)))
}

## ---------------------------------------------------------------------------
## Certified K = 11 strict-minimax constants (via package accessors)
## ---------------------------------------------------------------------------

mix_w <- LogisNormIntegral::minimax_weights(11L)
mix_s <- LogisNormIntegral::minimax_scales(11L)
stopifnot(abs(sum(mix_w) - 1) < 1e-12, all(mix_w > 0), all(mix_s > 0))

## Inner integrals, both vectorised over (eta, v) vectors of equal length.
inner_mixture <- function(eta, v) {
    z <- outer(eta, mix_s, "*") / sqrt(1 + outer(v, mix_s^2, "*"))
    as.numeric(pnorm(z) %*% mix_w)
}
inner_vvw <- function(eta, v) {
    p    <- plogis(eta)
    lam2 <- p * (1 - p) * (1 - 2 * p)
    p + 0.5 * v * lam2
}

## ---------------------------------------------------------------------------
## Data-generating process (identical to Example1.R)
## ---------------------------------------------------------------------------

common_pars <- list(
    theta_T = 0.5, theta_M = 1.5, theta_TM = 0.5, theta_C = -0.2,
    gamma_0 = 0.5, gamma_T = 1.2, gamma_C = 0.3, sigma_M = 0.8
)
n         <- 500L
N_sim     <- 2000L
R_qb      <- 4000L
scenarios <- list(
    list(name = "low",      theta_0 = -3.0),
    list(name = "moderate", theta_0 = -1.0),
    list(name = "high",     theta_0 =  0.0)
)

## Population NIE at the true coefficients: contrast z = 1, z' : 0 -> 1, c = 0.
## High-accuracy adaptive Gauss--Kronrod (`integrate`), the fixed estimand a
## frequentist interval must cover.
nie_population <- function(theta_0) {
    b      <- common_pars[["theta_M"]] + common_pars[["theta_TM"]] * 1
    a1     <- theta_0 + common_pars[["theta_T"]] * 1
    mu_M1  <- common_pars[["gamma_0"]] + common_pars[["gamma_T"]] * 1
    mu_M0  <- common_pars[["gamma_0"]]
    v      <- b^2 * common_pars[["sigma_M"]]^2
    I_ref  <- function(eta) {
        integrate(function(x) plogis(eta + sqrt(v) * x) * dnorm(x),
                  -Inf, Inf, rel.tol = 1e-12)[["value"]]
    }
    I_ref(a1 + b * mu_M1) - I_ref(a1 + b * mu_M0)
}

## One synthetic data set, fitted by ML. Returns the coefficient bundle the QB
## step needs, or NULL if the logistic fit did not converge.
simulate_and_fit <- function(theta_0, seed) {
    set.seed(seed)
    treat    <- rbinom(n, 1, 0.5)
    covar    <- rnorm(n)
    mediator <- common_pars[["gamma_0"]] + common_pars[["gamma_T"]] * treat +
        common_pars[["gamma_C"]] * covar +
        rnorm(n, sd = common_pars[["sigma_M"]])
    logit_p  <- theta_0 + common_pars[["theta_T"]] * treat +
        common_pars[["theta_M"]] * mediator +
        common_pars[["theta_TM"]] * treat * mediator +
        common_pars[["theta_C"]] * covar
    y  <- rbinom(n, 1, plogis(logit_p))
    df <- data.frame(y = y, treat = treat, mediator = mediator, covar = covar)

    mod_m <- lm(mediator ~ treat + covar, data = df)
    mod_y <- suppressWarnings(
        glm(y ~ treat * mediator + covar, data = df, family = binomial)
    )
    if (!isTRUE(mod_y[["converged"]])) return(NULL)

    list(
        theta_hat   = coef(mod_y),
        V_theta     = vcov(mod_y),
        gamma_hat   = coef(mod_m),
        V_gamma     = vcov(mod_m),
        sigma_M_hat = summary(mod_m)[["sigma"]],
        df_resid    = mod_m[["df.residual"]],
        py1         = mean(y)
    )
}

## 95% QB intervals from one fit, for both inner integrals off the same outer
## parameter draws (so the comparison isolates the inner-integral effect).
qb_intervals <- function(fit, seed) {
    set.seed(seed)
    theta_draws <- rmvnorm_chol(R_qb, fit[["theta_hat"]], fit[["V_theta"]])
    gamma_draws <- rmvnorm_chol(R_qb, fit[["gamma_hat"]], fit[["V_gamma"]])
    sigma_sq    <- fit[["df_resid"]] * fit[["sigma_M_hat"]]^2 /
        rchisq(R_qb, fit[["df_resid"]])

    z         <- 1
    covar_val <- 0
    mu_M1   <- as.numeric(gamma_draws %*% c(1, 1, covar_val))
    mu_M0   <- as.numeric(gamma_draws %*% c(1, 0, covar_val))
    alpha_1 <- as.numeric(theta_draws %*% c(1, z, 0, covar_val, 0))
    beta_M  <- theta_draws[, "mediator"] + z * theta_draws[, "treat:mediator"]
    eta_11  <- alpha_1 + beta_M * mu_M1
    eta_10  <- alpha_1 + beta_M * mu_M0
    v_both  <- beta_M^2 * sigma_sq

    ci_of <- function(inner) {
        yy  <- inner(c(eta_11, eta_10), c(v_both, v_both))
        nie <- yy[seq_len(R_qb)] - yy[R_qb + seq_len(R_qb)]
        unname(quantile(nie, c(0.025, 0.975)))
    }
    list(mixture = ci_of(inner_mixture), vvw = ci_of(inner_vvw))
}

covers <- function(ci, target) target >= ci[1] && target <= ci[2]

## ---------------------------------------------------------------------------
## Coverage sweep
## ---------------------------------------------------------------------------

set.seed(6489L)
data_seeds <- sample.int(.Machine[["integer.max"]], N_sim)
qb_seeds   <- sample.int(.Machine[["integer.max"]], N_sim)

run_scenario <- function(scn) {
    target <- nie_population(scn[["theta_0"]])
    cat(sprintf("--- Scenario: %-8s (theta_0 = %+.1f)  population NIE = %+.5f ---\n",
                scn[["name"]], scn[["theta_0"]], target))

    one_rep <- function(ds, qs) {
        fit <- simulate_and_fit(scn[["theta_0"]], ds)
        if (is.null(fit)) return(NULL)
        ci <- qb_intervals(fit, qs)
        list(
            py1     = fit[["py1"]],
            cov_mix = covers(ci[["mixture"]], target),
            cov_vvw = covers(ci[["vvw"]], target),
            w_mix   = ci[["mixture"]][2] - ci[["mixture"]][1],
            w_vvw   = ci[["vvw"]][2] - ci[["vvw"]][1]
        )
    }

    reps  <- Map(one_rep, data_seeds, qb_seeds)
    ok    <- Filter(Negate(is.null), reps)
    n_ok  <- length(ok)
    pick  <- function(field) vapply(ok, function(r) r[[field]], numeric(1L))

    cov_mix <- mean(pick("cov_mix"))
    cov_vvw <- mean(pick("cov_vvw"))
    se      <- function(p) sqrt(p * (1 - p) / n_ok)

    cat(sprintf("    replicates used: %d / %d   mean P(Y=1) = %.3f\n",
                n_ok, N_sim, mean(pick("py1"))))
    cat(sprintf("    minimax K=11 : coverage = %5.1f%% (SE %.1f)   mean width = %.4f\n",
                100 * cov_mix, 100 * se(cov_mix), mean(pick("w_mix"))))
    cat(sprintf("    V&VW Taylor  : coverage = %5.1f%% (SE %.1f)   mean width = %.4f\n\n",
                100 * cov_vvw, 100 * se(cov_vvw), mean(pick("w_vvw"))))

    data.frame(
        scenario   = scn[["name"]],
        theta_0    = scn[["theta_0"]],
        py1        = mean(pick("py1")),
        nie_pop    = target,
        n_used     = n_ok,
        cov_mix    = cov_mix,
        se_mix     = se(cov_mix),
        cov_vvw    = cov_vvw,
        se_vvw     = se(cov_vvw),
        width_mix  = mean(pick("w_mix")),
        width_vvw  = mean(pick("w_vvw"))
    )
}

cat("Repeated-sampling coverage of nominal 95% quasi-Bayesian NIE intervals\n")
cat("=========================================================================\n")
cat(sprintf("  n = %d per data set, N_sim = %d data sets, R_qb = %d outer draws.\n",
            n, N_sim, R_qb))
cat("  Target = fixed population NIE (true coefficients), contrast z = 1,\n")
cat("  z' : 0 -> 1, covariate c = 0. Inner integral: certified minimax K = 11\n")
cat("  vs Valeri--VanderWeele 2nd-order Taylor, off identical outer draws.\n\n")

summary_df <- Map(run_scenario, scenarios) |> Reduce(rbind, x = _)

cat("=========================================================================\n")
cat("Summary table\n\n")
print(summary_df, digits = 4L, row.names = FALSE)

cat("\nNotes:\n")
cat("  - Coverage is the fraction of the N_sim intervals containing the fixed\n")
cat("    population NIE. SE is the binomial Monte Carlo error sqrt(p(1-p)/N).\n")
cat("  - The minimax QB interval propagates full parameter uncertainty through\n")
cat("    the exact nonlinear map, so coverage tracks the nominal 95% level.\n")
cat("  - V&VW's O(v^2) Taylor bias shifts its intervals off the target. The\n")
cat("    resulting mis-coverage flips direction with outcome prevalence (severe\n")
cat("    under-coverage at low, over-coverage at moderate, and under-coverage at\n")
cat("    high), matching the sign-flipping bias in Example1.R and the manuscript\n")
cat("    method-comparison table.\n")

invisible(summary_df)
