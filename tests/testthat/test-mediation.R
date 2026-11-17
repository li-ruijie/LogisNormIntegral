# predict_indirect_effect_qb: single-profile contract and value ranges.

make_models <- function() {
    set.seed(42)
    N <- 200
    treat <- rbinom(N, 1, 0.5)
    covar <- rnorm(N)
    med <- 0.5 + 1.2 * treat + 0.3 * covar + rnorm(N, sd = 0.8)
    y <- rbinom(N, 1, plogis(-1 + 0.5 * treat + 1.5 * med - 0.2 * covar))
    df <- data.frame(y = y, treat = treat, med = med, covar = covar)
    list(mod_m = lm(med ~ treat + covar, data = df),
         mod_y = glm(y ~ treat + med + covar, data = df, family = binomial))
}

test_that("predict_indirect_effect_qb rejects multi-row covar_vals", {
    m <- make_models()
    expect_error(
        predict_indirect_effect_qb(m[["mod_y"]], m[["mod_m"]],
                                   covar_vals = data.frame(covar = c(0, 1)),
                                   n_sim = 50L),
        "single covariate profile"
    )
})

test_that("predict_indirect_effect_qb runs on a single covariate profile", {
    m <- make_models()
    set.seed(1)
    res <- predict_indirect_effect_qb(m[["mod_y"]], m[["mod_m"]],
                                      covar_vals = data.frame(covar = 0),
                                      n_sim = 400L)
    expect_true(is.finite(res[["NIE"]][["estimate"]]))
    expect_lt(res[["NIE"]][["ci_lower"]], res[["NIE"]][["ci_upper"]])
})
