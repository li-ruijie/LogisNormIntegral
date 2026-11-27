#!/usr/bin/env Rscript
## RandomIntercept.R
## ---------------------------------------------------------------------------
## Population-averaged probabilities of random-intercept logistic models on
## five public clustered data sets, by the certified mixture, by fixed-node
## Gauss-Hermite quadrature, and by the one-term probit scaling.
##
## With logit P(Y = 1 | x, u) = x'beta + u and u ~ N(0, tau^2), the probability
## of the outcome at covariates x, averaged over the clusters, is the
## logistic-normal integral I(x'beta, tau^2). Its variance argument v is the
## variance of the random intercept, which fitted models put well above the
## values met in causal mediation.
##
##   toenail        HSAUR3   onycholysis over seven visits, two treatments
##   respiratory    HSAUR3   respiratory status over four visits, two centres
##   ohio           geepack  wheeze at ages 7 to 10, maternal smoking
##   bacteria       MASS     H. influenzae carriage in children with otitis media
##   Contraception  mlmRev   contraceptive use by district, Bangladesh
##
## Each model is fitted by lme4::glmer with adaptive Gauss-Hermite quadrature
## at 50 points. The evaluators are then compared on the fixed-effect linear
## predictors of the fit, with adaptive Gauss-Kronrod quadrature as the
## reference, itself checked against a 500-node Gauss-Hermite rule.
##
## Two one-term approximations are included, the attenuation formula of
## Zeger, Liang, and Albert (1988), plogis(eta / sqrt(1 + c^2 v)) with
## c = 16 sqrt(3) / (15 pi), and the probit form pnorm(eta / sqrt(1.6^2 + v)).
##
## Usage: Rscript RandomIntercept.R
## Requires lme4 and statmod, and HSAUR3, geepack, MASS, and mlmRev for their
## data alone.
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
needed  <- c("lme4", "statmod", "HSAUR3", "geepack", "MASS", "mlmRev")
missing <- Filter(function(p) !requireNamespace(p, quietly = TRUE), needed)
if (length(missing) > 0L) {
    stop("Install the packages: ", paste(missing, collapse = ", "))
}

## A data set of a package that need not be attached.
package_data <- function(name, package) {
    env <- new.env()
    utils::data(list = name, package = package, envir = env)
    as.data.frame(env[[name]])
}
set_col <- function(d, name, value) `[[<-`(d, name, value = value)

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
    Map(gh_eval, c(GH_10 = 10L, GH_20 = 20L, GH_30 = 30L, GH_50 = 50L, GH_100 = 100L)),
    Map(mixture, c(mixture_K11 = 11L, mixture_K22 = 22L))
)
gh_500 <- gh_eval(500L)

fit_model <- function(formula, d) {
    lme4::glmer(formula, data = d, family = binomial, nAGQ = 50L,
                control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
}
intercept_variance <- function(fit) sum(as.data.frame(lme4::VarCorr(fit))[["vcov"]])

## Compare the evaluators on the distinct fixed-effect linear predictors.
compare <- function(label, fit) {
    v     <- intercept_variance(fit)
    eta   <- unique(round(as.vector(model.matrix(fit) %*% lme4::fixef(fit)), 10))
    truth <- reference(eta, v)
    data.frame(
        data      = label,
        n         = nrow(model.frame(fit)),
        clusters  = unname(lme4::ngrps(fit)),
        v         = v,
        eta_min   = min(eta),
        eta_max   = max(eta),
        reference = max(abs(truth - gh_500(eta, v)))
    ) |> cbind(as.list(mapply(function(f) max(abs(f(eta, v) - truth)), evaluators)))
}

toenail     <- package_data("toenail", "HSAUR3")
toenail     <- set_col(toenail, "y", as.numeric(toenail[["outcome"]] == "moderate or severe"))
respiratory <- package_data("respiratory", "HSAUR3")
respiratory <- respiratory[respiratory[["month"]] != "0", ]
respiratory <- respiratory |>
    set_col("y", as.numeric(respiratory[["status"]] == "good")) |>
    set_col("age10", respiratory[["age"]] / 10)
ohio        <- package_data("ohio", "geepack")
bacteria    <- package_data("bacteria", "MASS")
bacteria    <- bacteria |>
    set_col("y", as.numeric(bacteria[["y"]] == "y")) |>
    set_col("late", as.numeric(bacteria[["week"]] > 2))
contra      <- package_data("Contraception", "mlmRev")
contra      <- contra |>
    set_col("y", as.numeric(contra[["use"]] == "Y")) |>
    set_col("age10", contra[["age"]] / 10)

fits <- list(
    toenail       = fit_model(y ~ treatment * time + (1 | patientID), toenail),
    respiratory   = fit_model(y ~ centre + treatment + gender + age10 + (1 | subject), respiratory),
    ohio          = fit_model(resp ~ age + smoke + (1 | id), ohio),
    bacteria      = fit_model(y ~ trt + late + (1 | ID), bacteria),
    Contraception = fit_model(y ~ age10 + I(age10^2) + urban + livch + (1 | district), contra)
)
results <- Map(compare, names(fits), fits) |> Reduce(f = rbind)

cat("Population-averaged probabilities of random-intercept logistic models\n")
cat("=========================================================================\n")
cat("  v is the fitted variance of the random intercept. Each evaluator column is\n")
cat("  the largest absolute error over the fixed-effect linear predictors of the\n")
cat("  fit, against adaptive quadrature. The column 'reference' is the largest\n")
cat("  difference between that quadrature and a 500-node Gauss-Hermite rule.\n\n")
options(width = 200)
## Errors to two significant digits, rounded here from full precision.
error_cols <- c("reference", names(evaluators))
shown      <- Reduce(function(d, name) set_col(d, name, signif(d[[name]], 2L)), error_cols, results)
print(format(shown, digits = 5L), row.names = FALSE)

## The contrast an analyst reports from the toenail trial: the difference
## between the arms in the population-averaged probability of the outcome.
fit      <- fits[["toenail"]]
v        <- intercept_variance(fit)
beta     <- lme4::fixef(fit)
months   <- c(0, 3, 6, 12)
eta_arm  <- function(second) {
    beta[[1L]] + second * beta[[2L]] + months * (beta[["time"]] + second * beta[[4L]])
}
contrast <- function(f) f(eta_arm(1), v) - f(eta_arm(0), v)
table    <- Map(contrast, c(list(reference = reference), evaluators)) |> Reduce(f = rbind)
table    <- `dimnames<-`(table, list(c("reference", names(evaluators)), paste0("month_", months)))
cat(sprintf("\nToenail trial (v = %.3f): %s minus %s, the difference in the\n", v,
            levels(toenail[["treatment"]])[[2L]], levels(toenail[["treatment"]])[[1L]]))
cat("population-averaged probability of moderate or severe onycholysis:\n\n")
print(table, digits = 8L)
cat("\nPopulation-averaged probabilities in the two arms, by the reference:\n\n")
arms <- rbind(reference(eta_arm(0), v), reference(eta_arm(1), v))
print(`dimnames<-`(arms, list(levels(toenail[["treatment"]]), paste0("month_", months))), digits = 5L)
cat("\nRelative error of each contrast, in per cent:\n\n")
print(round(100 * sweep(table[-1L, , drop = FALSE], 2L, table[1L, ], function(x, r) x / r - 1), 3L))

invisible(results)
