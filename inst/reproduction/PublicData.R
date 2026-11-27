#!/usr/bin/env Rscript
## PublicData.R
## ---------------------------------------------------------------------------
## The odds ratio of the natural indirect effect on a public data set, exact by
## the certified mixture against the rare-outcome closed form of Valeri and
## VanderWeele (2013).
##
## Example1.R makes the same comparison at synthetic coefficients, where one
## residual standard deviation of the mediator multiplies the odds of the
## outcome by about five (v = 2.56). This script measures how far the two
## answers lie apart at the mediator strength of a fitted model.
##
##   WCGS       Western Collaborative Group Study (Rosenman et al. 1975), from
##              the CRAN package faraway. Exposure: body-mass index of 25 or
##              more. Mediator: systolic blood pressure. Outcome: coronary heart
##              disease.
##
## The analysis fits a normal linear model for the mediator and a logistic
## model for the outcome with an exposure-by-mediator interaction, as in
## Section 5 of the accompanying paper. Body-mass index and blood pressure are
## states rather than interventions, and both were measured at the same visit,
## so the computation evaluates the mediation formula at the fitted models and
## offers no causal reading of the data.
##
## The closed form is one number for the data set, exp{(theta_M + theta_TM)
## gamma_T}. The exact odds ratio is conditional on the covariates and varies
## with them, so its geometric mean, minimum, and maximum over the individuals
## are reported.
##
## Usage: Rscript PublicData.R
## Requires the package faraway for its data alone.
##
## R style: 4-space indentation, Map/Reduce/Filter instead of for-loops,
## `[[` instead of `$`, base pipe `|>`, functional replacement forms.

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

## A data set of a package that need not be attached.
package_data <- function(name, package) {
    if (!requireNamespace(package, quietly = TRUE)) {
        stop(sprintf("Install the package %s, which holds the %s data.", package, name))
    }
    env <- new.env()
    utils::data(list = name, package = package, envir = env)
    as.data.frame(env[[name]])
}
set_col <- function(d, name, value) `[[<-`(d, name, value = value)
odds <- function(p) p / (1 - p)

## Adaptive quadrature, the independent check on the mixture.
reference_integral <- function(eta, v) {
    integrate(function(m) plogis(eta + sqrt(v) * m) * dnorm(m), -Inf, Inf,
              rel.tol = 1e-12, abs.tol = 1e-14)[["value"]]
}

## Fit the two models and compare the two odds ratios of the indirect effect
## of moving the exposure from 0 to 1 in the mediator, with the direct exposure
## held at 1.
compare <- function(label, d, a, m, y, covars) {
    d       <- d[complete.cases(d[, c(a, m, y, covars)]), ]
    rhs     <- paste(covars, collapse = " + ")
    fit_m   <- lm(as.formula(paste(m, "~", a, "+", rhs)), data = d)
    fit_y   <- glm(as.formula(paste(y, "~", a, "*", m, "+", rhs)), data = d, family = binomial)
    gamma   <- coef(fit_m)
    theta   <- coef(fit_y)
    sigma_M <- summary(fit_m)[["sigma"]]
    inter   <- grep(":", names(theta), value = TRUE)
    beta_M  <- function(z) theta[[m]] + theta[[inter]] * z
    ## Linear predictors of each individual with the exposure terms removed.
    mean_m0 <- as.vector(predict(fit_m, newdata = set_col(d, a, 0)))
    x_y     <- model.matrix(fit_y)
    rest    <- !(colnames(x_y) %in% c(a, m, inter))
    base_y  <- as.vector(x_y[, rest, drop = FALSE] %*% theta[rest])
    eta     <- function(z_m) base_y + theta[[a]] + beta_M(1) * (mean_m0 + gamma[[a]] * z_m)
    v       <- beta_M(1)^2 * sigma_M^2
    mixture <- function(z_m) LogisNormIntegral::logistic_normal_integral(eta(z_m), v, K = 22L)
    quad    <- function(z_m) mapply(function(e) reference_integral(e, v), eta(z_m))
    p_11    <- mixture(1)
    p_10    <- mixture(0)
    or_i    <- odds(p_11) / odds(p_10)
    data.frame(
        data          = label,
        n             = nrow(d),
        prevalence    = mean(d[[y]]),
        v_exposed     = v,
        v_unexposed   = beta_M(0)^2 * sigma_M^2,
        OR_rare_form  = exp(beta_M(1) * gamma[[a]]),
        OR_exact_mean = exp(mean(log(or_i))),
        OR_exact_min  = min(or_i),
        OR_exact_max  = max(or_i),
        rare_over_exact = exp(beta_M(1) * gamma[[a]]) / exp(mean(log(or_i))),
        mixture_abs_err = max(abs(p_11 - quad(1)), abs(p_10 - quad(0)))
    )
}

wcgs <- package_data("wcgs", "faraway")
wcgs <- wcgs |>
    set_col("overweight", as.numeric(703 * wcgs[["weight"]] / wcgs[["height"]]^2 >= 25)) |>
    set_col("type_a", as.numeric(wcgs[["dibep"]] == "A")) |>
    set_col("chd01", as.numeric(wcgs[["chd"]] == "yes"))

results <- compare("WCGS", wcgs, "overweight", "sdp", "chd01", c("age", "cigs", "type_a"))

cat("The odds ratio of the natural indirect effect on a public data set\n")
cat("=========================================================================\n")
cat("  Exact: the certified K = 22 mixture, conditional on the covariates of each\n")
cat("  individual. Closed form: Valeri and VanderWeele (2013), rare outcome.\n\n")
options(width = 200)
print(results, digits = 5L, row.names = FALSE)
cat("\nNotes:\n")
cat("  - v is the squared log odds ratio of the outcome for one residual standard\n")
cat("    deviation of the mediator, (theta_M + theta_TM z)^2 sigma_M^2.\n")
cat("  - mixture_abs_err is the largest difference between the mixture and adaptive\n")
cat("    quadrature over every potential-outcome mean evaluated.\n")

invisible(results)
