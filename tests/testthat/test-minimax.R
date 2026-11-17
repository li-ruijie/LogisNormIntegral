# Certified minimax mixture: structure, the uniform sup-norm bound, the
# v = 0 reduction, and the default K = 22 signed frontier.

positive_K <- c(1:11, 14L, 15L, 16L, 17L)   # positive-class certs (all weights > 0)
z_grid <- seq(-8, 8, by = 0.25)

test_that("positive-class certs are well-formed", {
    invisible(Map(
        function(K) {
            w <- minimax_weights(K)
            s <- minimax_scales(K)
            expect_length(w, K)
            expect_length(s, K)
            expect_equal(sum(w), 1, tolerance = 1e-9)
            expect_true(all(w > 0))
            expect_true(all(s > 0))
            expect_equal(length(unique(s)), K)
        },
        positive_K
    ))
})

test_that("mixture stays within its certified sup-norm bound on a grid", {
    # Grid maximum can only under-estimate the continuous sup, which is the
    # certified minimax, so |mixture - Lambda| must not exceed Delta_K. The pad
    # absorbs double-precision evaluation noise, measured at up to 4.5e-16 on a
    # 60001-point grid, and is absolute because that noise is. A relative pad
    # would shrink with Delta_K and fall below the noise at the small-Delta end.
    invisible(Map(
        function(K) {
            err <- max(abs(phi_minimax_mixture(z_grid, K) - plogis(z_grid)))
            expect_lte(err, minimax_error(K) + 1e-14)
        },
        positive_K
    ))
})

test_that("logistic_normal_integral at v = 0 reduces to the mixture CDF", {
    eta <- c(-2, -0.5, 0, 0.5, 2)
    invisible(Map(
        function(K) expect_equal(logistic_normal_integral(eta, 0, K),
                                 phi_minimax_mixture(eta, K)),
        c(1L, 8L, 11L, 22L)
    ))
})

test_that("K = 22 signed frontier is usable and below 1e-13", {
    expect_true(minimax_constants[["K22"]][["signed"]])
    expect_equal(sum(minimax_weights(22)), 1, tolerance = 1e-9)
    err <- max(abs(phi_minimax_mixture(z_grid, 22) - plogis(z_grid)))
    # The threshold is permissive. The grid maximum is about 4e-15, the 3.6e-15
    # mathematical bound plus fp64 evaluation noise.
    expect_lte(err, 1e-13)
})

test_that("certified error constants match the proved values", {
    expect_equal(minimax_error(11), 1.2584376967160870e-11, tolerance = 1e-12)
    expect_equal(minimax_error(8), 2.1086e-9, tolerance = 1e-4)
})

test_that("invalid K is rejected", {
    expect_error(minimax_weights(12))
    expect_error(phi_minimax_mixture(0, K = 0))
})

test_that("phi_small_v order-2 correction strictly improves on order-1", {
    # Reference matching phi_small_v's convention: phi_small_v(eta, 0, 0) =
    # 1 / (1 + e^eta) = plogis(-eta), so the integrand is plogis(-x) on N(eta, v).
    I_true <- function(eta, v) {
        integrate(function(x) plogis(-x) * dnorm(x, eta, sqrt(v)),
                  -Inf, Inf, rel.tol = 1e-13, abs.tol = 1e-15)[["value"]]
    }
    eta <- 0.7
    invisible(Map(
        function(v) {
            truth <- I_true(eta, v)
            e0 <- abs(phi_small_v(eta, v, 0) - truth)
            e1 <- abs(phi_small_v(eta, v, 1) - truth)
            e2 <- abs(phi_small_v(eta, v, 2) - truth)
            # Each asymptotic order strictly improves in the small-v regime.
            # The pre-fix finite-difference order-2 term failed this (e2 > e1).
            expect_lt(e1, e0)
            expect_lt(e2, e1)
        },
        c(0.125, 0.25, 0.5)
    ))
    # Closed-form 4th derivative gives ~v^3 accuracy, while the old h^4 stencil
    # gave ~8e-3 here (noise), so this bound is the regression guard.
    expect_lt(abs(phi_small_v(eta, 0.125, 2) - I_true(eta, 0.125)), 1e-5)
})

test_that("the certified accessors work on a loaded-but-unattached namespace", {
    ## Regression: minimax_constants and minimax_delta live in the LazyData
    ## environment, which library() attaches to the search path but which is not
    ## on a package function's namespace lookup chain. Bare references therefore
    ## resolved only when the caller happened to attach the package, so every
    ## exported function touching the certified constants failed under `::`,
    ## including from the package's own reproduction scripts.
    skip_on_cran()
    script <- tempfile(fileext = ".R")
    writeLines(c(
        'v <- c(LogisNormIntegral::minimax_error(11L),',
        '       LogisNormIntegral::minimax_weights(11L)[1],',
        '       LogisNormIntegral::minimax_scales(11L)[1],',
        '       LogisNormIntegral::logistic_normal_integral(0, 1),',
        '       LogisNormIntegral::phi_minimax_mixture(0))',
        'stopifnot(all(is.finite(v)))',
        'cat("UNATTACHED_OK\n")'
    ), script)
    out <- suppressWarnings(system2(
        file.path(R.home("bin"), "Rscript"),
        c("--vanilla", shQuote(script)), stdout = TRUE, stderr = TRUE))
    expect_true(any(grepl("UNATTACHED_OK", out)),
                info = paste(out, collapse = "\n"))
})
