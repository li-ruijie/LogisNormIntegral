# Asymptotic series: the scalar contract for the dispatch functions and the
# vectorised phi_small_v. Accuracy of phi_small_v(order = 2) is covered in
# test-minimax.R. The large-v series is asymptotic only (it diverges near the
# domain boundaries and for moderate v), so its accuracy is not asserted here.

test_that("scalar asymptotic functions reject vector input with a clear error", {
    expect_error(phi_large_v_base(c(3, 4), 10), "scalar")
    expect_error(phi_large_v_global(c(1, 2), 8), "scalar")
    expect_error(phi_analytic(c(1, 2), c(2, 8)), "scalar")
})

test_that("scalar asymptotic functions still evaluate on scalar input", {
    expect_length(phi_large_v_base(3, 12), 1)
    expect_length(phi_large_v_global(2, 30), 1)
    expect_length(phi_analytic(0.5, 2), 1)
})

test_that("phi_small_v is vectorised over eta and v", {
    eta <- c(-1, 0, 1)
    v <- c(0.5, 1, 2)
    out <- phi_small_v(eta, v, order = 2)
    expect_length(out, 3)
    expect_equal(out, mapply(function(e, vv) phi_small_v(e, vv, order = 2), eta, v))
})

test_that("phi_analytic selects a branch and never silently returns a non-probability", {
    ## the cell that motivated the adaptive rule: the fixed v <= 6 threshold
    ## handed this to the large-v series, which returned about 114.
    v <- suppressWarnings(phi_analytic(1, 6.5, order = 2))
    expect_gte(as.numeric(v), 0)
    expect_lte(as.numeric(v), 1)
    expect_true(attr(v, "branch") %in% c("small", "large"))
    expect_warning(phi_analytic(1, 6.5, order = 2), "asymptotic series is not reliable")
    expect_silent(phi_analytic(0.5, 1, order = 2))
    ## Regression: the unit-interval rule can reject a branch only in favour of
    ## an in-range one, so it cannot fire when both branches have diverged, and
    ## the ratio test is silent there. This cell returned -0.226 against a true
    ## 0.212 with no warning at all, the worst of the three such arguments on the
    ## 119400-point measurement grid behind the documented guarantee. The range
    ## test on the returned value is what makes this test's title true.
    expect_warning(phi_analytic(3, 11, order = 2), "left the unit interval")
})

test_that("the warning is judged on the branch returned, not on the rejected one", {
    ## Regression: the warning tested min(r_small, r_large), which is the
    ## diagnostic of the returned branch only while that branch is the argmin.
    ## The unit-interval rule breaks exactly that. At (3.5, 12) order 2 the
    ## small branch returns -0.289 and is rejected at ratio 0.056, so the large
    ## branch is returned at ratio 4.66 and the minimum silenced the warning.
    ## The value was 0.719 against a true 0.184, an error of 0.53 on a quantity
    ## confined to [0, 1].
    expect_warning(phi_analytic(3.5, 12, order = 2), "asymptotic series is not reliable")
    expect_warning(phi_analytic(7, 9, order = 2), "asymptotic series is not reliable")
    val <- suppressWarnings(phi_analytic(3.5, 12, order = 2))
    expect_identical(attr(val, "branch"), "large")
})

test_that("forcing a branch reproduces the underlying series exactly", {
    expect_equal(as.numeric(phi_analytic(0.5, 1, 2, method = "small")), phi_small_v(0.5, 1, 2))
    expect_equal(as.numeric(phi_analytic(10, 20, 2, method = "large")),
                 phi_large_v_global(10, 20, 2))
})

test_that("a forced branch runs the same warn logic as the adaptive path", {
    ## Regression: `method = "small"` and `method = "large"` returned before the
    ## warn block, so the documented `warn` contract was inoperative on both. At
    ## (3.5, 12) order 2 the forced large branch returns 0.719 against a true
    ## 0.184 at a first-omitted-term ratio of 4.66, which the adaptive path had
    ## warned on since the selected-branch correction.
    expect_warning(phi_analytic(3.5, 12, order = 2, method = "large"),
                   "asymptotic series is not reliable")
    expect_warning(phi_analytic(3.5, 12, order = 2, method = "large"), "large-v branch")
    expect_silent(phi_analytic(3.5, 12, order = 2, method = "large", warn = FALSE))
    ## Regression: the ratio test alone cannot catch a diverged truncation. The
    ## same forced small branch returns -0.289 at a ratio of 0.056, below the
    ## threshold, so only the unit-interval test raises this one.
    expect_warning(phi_analytic(3.5, 12, order = 2, method = "small"),
                   "left the unit interval")
    ## A forced branch whose truncation is still well behaved stays silent.
    expect_silent(phi_analytic(0.5, 1, order = 2, method = "small"))
})

test_that("the asymptotic functions reject an order they do not implement", {
    ## Regression: order = 3 read past the end of the four-entry .small_terms
    ## and .large_terms tables, so .trunc_ratio returned NA and phi_analytic
    ## died on "missing value where TRUE/FALSE needed". The series themselves
    ## stop at the `order >= 2` branch, so any higher order also returned the
    ## order-2 value with no signal.
    expect_error(phi_analytic(0.5, 1, order = 3), "order")
    expect_error(phi_small_v(0.5, 1, order = 3), "order")
    expect_error(phi_large_v_base(3, 12, order = 3), "order")
    expect_error(phi_large_v_global(2, 30, order = 3), "order")
    expect_error(phi_analytic(0.5, 1, order = -1), "order")
    expect_error(phi_analytic(0.5, 1, order = 1.5), "order")
    expect_error(phi_analytic(0.5, 1, order = c(1, 2)), "order")
})

test_that("every implemented order still evaluates", {
    invisible(Map(
        function(o) {
            expect_length(suppressWarnings(phi_analytic(0.5, 1, order = o)), 1)
            expect_length(phi_small_v(0.5, 1, order = o), 1)
            expect_length(phi_large_v_global(2, 30, order = o), 1)
        },
        0:2
    ))
})
