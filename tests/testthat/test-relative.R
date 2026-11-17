# The opt-in relative-accuracy evaluation: the exponential-tilt identity
# I(eta, v) = exp(eta + v / 2) * I(-eta - v, v), applied for eta < -v / 2, and
# the relative bound Delta_K / I(-v / 2, v) that it gives at every eta.

# Adaptive quadrature split at the midpoint of the logistic factor, so that the
# exponential tail is resolved to relative and not only absolute accuracy.
quad_reference <- function(eta, v) {
    s <- sqrt(v)
    f <- function(z) plogis(eta + s * z) * dnorm(z)
    mid <- min(max(-eta / s, -37), 37)
    piece <- function(a, b) {
        integrate(f, a, b, rel.tol = 1e-13, abs.tol = 0, subdivisions = 4000L)[["value"]]
    }
    piece(-38, mid) + piece(mid, 38)
}

rel_grid <- expand.grid(eta = seq(-60, 30, by = 5), v = c(1, 4, 9, 16))

test_that("the default evaluation is the plain mixture sum, unchanged", {
    plain <- function(eta, v, K) {
        w <- minimax_weights(K)
        s <- minimax_scales(K)
        Reduce(`+`, Map(function(wi, si) wi * pnorm(si * eta / sqrt(1 + si^2 * v)), w, s), 0)
    }
    invisible(Map(
        function(K) {
            expect_identical(
                logistic_normal_integral(rel_grid[["eta"]], rel_grid[["v"]], K),
                plain(rel_grid[["eta"]], rel_grid[["v"]], K)
            )
            expect_identical(
                logistic_normal_integral(rel_grid[["eta"]], rel_grid[["v"]], K, relative = FALSE),
                logistic_normal_integral(rel_grid[["eta"]], rel_grid[["v"]], K)
            )
        },
        c(8L, 11L, 22L)
    ))
})

test_that("relative = TRUE leaves arguments at or above -v / 2 untouched", {
    keep <- rel_grid[["eta"]] >= -rel_grid[["v"]] / 2
    expect_identical(
        logistic_normal_integral(rel_grid[["eta"]][keep], rel_grid[["v"]][keep], 11, relative = TRUE),
        logistic_normal_integral(rel_grid[["eta"]][keep], rel_grid[["v"]][keep], 11)
    )
})

test_that("the tilted evaluation stays within the relative bound on a grid", {
    ref <- mapply(quad_reference, rel_grid[["eta"]], rel_grid[["v"]])
    invisible(Map(
        function(K) {
            tilted <- logistic_normal_integral(rel_grid[["eta"]], rel_grid[["v"]], K, relative = TRUE)
            bound <- minimax_relative_error(rel_grid[["v"]], K)
            # The pad covers the reference, which is accurate to about 1e-13 relative.
            expect_true(all(abs(tilted - ref) / ref <= bound + 1e-12))
        },
        c(8L, 11L)
    ))
})

test_that("the plain sum loses relative accuracy where the tilted one keeps it", {
    ref <- quad_reference(-30, 1)
    plain <- logistic_normal_integral(-30, 1, K = 11)
    tilted <- logistic_normal_integral(-30, 1, K = 11, relative = TRUE)
    expect_gt(abs(plain - ref) / ref, 1e-2)
    expect_lt(abs(tilted - ref) / ref, 1e-10)
})

test_that("minimax_relative_error is the certified ratio and is infinite without a guarantee", {
    v <- c(1, 9, 36)
    floor_value <- logistic_normal_integral(-v / 2, v, K = 11) - minimax_error(11)
    expect_equal(minimax_relative_error(v, K = 11), minimax_error(11) / floor_value)
    expect_lt(minimax_relative_error(1, K = 11), 3.2e-11)
    expect_lt(minimax_relative_error(36, K = 11), 6.1e-9)
    expect_identical(minimax_relative_error(2000, K = 11), Inf)
    expect_error(minimax_relative_error(1, K = 12), "K must be")
})

test_that("relative = TRUE recycles and validates like the default", {
    expect_length(logistic_normal_integral(c(-40, -20, 0), 4, K = 11, relative = TRUE), 3L)
    expect_length(logistic_normal_integral(-40, c(1, 4, 9), K = 11, relative = TRUE), 3L)
    expect_error(logistic_normal_integral(0, 1, K = 12, relative = TRUE), "K must be")
    expect_error(logistic_normal_integral(0, 1, K = 11, relative = NA), "relative")
})
