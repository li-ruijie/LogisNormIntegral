#' @title Closed-Form Asymptotic and Uniform Solutions for the Logistic-Normal Integral
#' @name logistic_normal
#'
#' @description
#' The logistic-normal integral lacks an exact closed-form solution
#' over the continuous real line. This package provides two classes of solutions:
#'
#' 1. **Mathematical Asymptotics (Exact Duality)**: Duality between the logistic-
#'    normal integral and Zwegers' Mordell integral yields exact asymptotic Taylor
#'    expansions. Accuracy varies across the parameter space (notably at the
#'    crossover variance v ~ 6.0), which can destabilise confidence-interval bounds.
#'
#' 2. **Uniform Minimax Normal Scale Mixture**: A K-term minimax mixture
#'    (Monahan & Stefanski 1992, certified in this project for the strict
#'    K = 1, ..., 11, the relaxed Bernstein-Akhiezer K = 14, ..., 17, and the
#'    signed K = 22 frontier)
#'    gives a uniform absolute error <= Delta_K* across the entire parameter space.
#'    For K = 1, ..., 11 these are the certified global minimax, proved to
#'    equioscillate at 2K points and to be the unique global minimiser via the
#'    varisolvency of the mixture family (see the package vignette). Monahan &
#'    Stefanski conjectured the equioscillation and computed the K <= 8 values
#'    numerically. The error bound Delta_K* is a verified property of the shipped
#'    mixture, independent of the optimality claim.
#'
#' @details
#' The certified minimax errors (Delta_K*) are:
#' \tabular{ll}{
#'   K = 1 \tab 9.46e-3 \cr
#'   K = 2 \tab 5.10e-4 \cr
#'   K = 3 \tab 4.36e-5 \cr
#'   K = 4 \tab 4.72e-6 \cr
#'   K = 5 \tab 5.95e-7 \cr
#'   K = 6 \tab 8.39e-8 \cr
#'   K = 7 \tab 1.29e-8 \cr
#'   K = 8 \tab 2.11e-9 \cr
#'   K = 9 \tab 3.66e-10 \cr
#'   K = 10 \tab 6.65e-11 \cr
#'   K = 11 \tab 1.26e-11 \cr
#'   K = 14 (BA-(K+1) cert) \tab 7.23e-12 \cr
#'   K = 15 (BA-(K+1) cert) \tab 7.21e-12 \cr
#'   K = 16 (BA-(K+1) cert) \tab 4.07e-12 \cr
#'   K = 17 (BA-K cert) \tab 2.54e-12 \cr
#'   K = 22 (signed frontier) \tab 3.60e-15
#' }
#' The K = 14, ..., 17 entries are relaxed Bernstein-Akhiezer certificates of
#' the positive class, verified bounds on the shipped mixtures rather than
#' proved global minimax values.
#' The K = 22 entry is the best certified bound, in the signed class
#' (\eqn{\sum_i w_i = 1} with one negative atom). See `minimax_constants$K22`.
#' This table is shown to three significant figures, while `minimax_error()` returns
#' the certified value at full double precision. Evaluating the mixture in
#' double-precision arithmetic adds rounding at the level of machine epsilon,
#' which sits far below every tabulated bound for \eqn{K \le 17}. At the signed
#' K = 22 frontier the accumulated evaluation noise is about one eighth of the
#' bound, so a double-precision evaluation can exceed it by a few times
#' \eqn{10^{-16}} while the exact error stays below it.
#' See `vignette("certified-minimax")` for the equioscillation, local-minimality,
#' and global-uniqueness (varisolvency) results.
#'
#' @return For `phi_minimax_mixture()` / `logistic_normal_integral()`: the
#'   approximated integral value. For accessors: the underlying certified
#'   weights, scales, or error bound.
NULL

#' Evaluate the certified K-term minimax mixture CDF approximation
#'
#' Computes \eqn{f^\star(z) = \sum_i w_i^\star \Phi(s_i^\star z)} using the
#' certified weights and scales for \eqn{K \in \{1, \ldots, 11\}} (strict), the relaxed
#' Bernstein-Akhiezer certs at K = 14, ..., 17, or the signed
#' frontier at K = 22 (the best certified bound, \eqn{\Delta \le 3.605 \times 10^{-15}}).
#'
#' @param z Numeric (scalar or vector) argument.
#' @param K Integer in 1..11 (strict), 14..17 (BA-cert), or 22 (signed
#'   frontier), giving the approximation order. Default 22 (signed frontier, the
#'   lowest certified error).
#' @return The mixture CDF evaluated at z.
#' @export
#' @examples
#' phi_minimax_mixture(0, K = 8)      # ~ 0.5
#' phi_minimax_mixture(2, K = 10)     # close to Lambda(2) = 1/(1+exp(-2))
phi_minimax_mixture <- function(z, K = 22L) {
  K <- as.integer(K)
  if (!(K %in% c(1:11, 14L, 15L, 16L, 17L, 22L))) {
    stop("K must be 1..11 (strict certified), 14..17 (BA-cert), or 22 (certified signed frontier).")
  }
  w <- minimax_weights(K)
  s <- minimax_scales(K)
  result <- numeric(length(z))
  for (i in seq_len(K)) {
    result <- result + w[i] * stats::pnorm(s[i] * z)
  }
  result
}

#' Evaluate the logistic-normal integral via the certified minimax mixture
#'
#' Computes
#' \deqn{I(\eta, v) = E[\Lambda(Y)], \quad Y \sim N(\eta, v)}
#' using the K-term certified minimax mixture approximation to \eqn{\Lambda}.
#'
#' The error is bounded by \eqn{|I - \widehat{I}_K| \le \Delta_K^\star},
#' where \eqn{\Delta_K^\star} is the certified bound returned by
#' `minimax_error()` (for K <= 11 the proved global minimax).
#'
#' @param eta Mean of the Gaussian predictor (scalar or vector).
#' @param v Variance of the Gaussian predictor (scalar or vector).
#' @param K Integer in 1..11 (strict), 14..17 (BA-cert), or 22 (signed
#'   frontier), giving the approximation order. Default 22 (signed frontier, the
#'   lowest certified error).
#' @return The approximated integral, vectorised over (eta, v).
#' @export
#' @examples
#' # E[Lambda(Y)], Y ~ N(eta, v), via the certified mixture
#' logistic_normal_integral(0, 1)                  # 0.5 by symmetry
#' logistic_normal_integral(c(-1, 0, 1), v = 2, K = 11)
logistic_normal_integral <- function(eta, v, K = 22L) {
  K <- as.integer(K)
  if (!(K %in% c(1:11, 14L, 15L, 16L, 17L, 22L))) {
    stop("K must be 1..11 (strict certified), 14..17 (BA-cert), or 22 (certified signed frontier).")
  }
  w <- minimax_weights(K)
  s <- minimax_scales(K)
  result <- numeric(length(eta))
  for (i in seq_len(K)) {
    scaled_eta <- s[i] * eta / sqrt(1 + s[i]^2 * v)
    result <- result + w[i] * stats::pnorm(scaled_eta)
  }
  result
}

#' Certified minimax error bound Delta_K*
#'
#' Returns the certified error bound \eqn{\Delta_K^\star} for the K-term mixture
#' approximation of the logistic CDF. Source: the audit deposit in `inst/certs/`.
#'
#' @param K Integer in 1..11 (strict), 14..17 (BA-cert), or 22 (signed frontier).
#' @return The certified \eqn{\Delta_K^\star} (for K <= 11 the proved global
#'   minimax).
#' @export
#' @examples
#' minimax_error(8)   # 2.1e-9
#' minimax_error(11)  # 1.26e-11
minimax_error <- function(K = 22L) {
  K <- as.integer(K)
  if (!(K %in% c(1:11, 14L, 15L, 16L, 17L, 22L))) {
    stop("K must be 1..11 (strict certified), 14..17 (BA-cert), or 22 (certified signed frontier).")
  }
  ## Qualified, not bare. The datasets live in the LazyData environment, which
  ## `library()` attaches to the search path but which is not on a package
  ## function's namespace-imports-base lookup chain. A bare reference therefore
  ## resolves only when the caller happened to attach the package, and fails
  ## under `LogisNormIntegral::minimax_error()` on a loaded-but-unattached
  ## namespace, which is how the reproduction scripts call it.
  LogisNormIntegral::minimax_delta[[paste0("K", K)]]
}

#' Certified minimax weights for K-term mixture
#'
#' @param K Integer in 1..11 (strict), 14..17 (BA-cert), or 22 (signed frontier).
#' @return Numeric vector of K weights, summing to 1, all positive for K <= 17.
#'   The signed K = 22 frontier carries one negative atom.
#' @export
#' @examples
#' w <- minimax_weights(8)
#' sum(w)   # 1
minimax_weights <- function(K = 22L) {
  K <- as.integer(K)
  if (!(K %in% c(1:11, 14L, 15L, 16L, 17L, 22L))) {
    stop("K must be 1..11 (strict certified), 14..17 (BA-cert), or 22 (certified signed frontier).")
  }
  ## Qualified rather than bare, for the LazyData lookup reason given in
  ## `minimax_error()`.
  LogisNormIntegral::minimax_constants[[paste0("K", K)]]$weights
}

#' Certified minimax scales for K-term mixture
#'
#' @param K Integer in 1..11 (strict), 14..17 (BA-cert), or 22 (signed frontier).
#' @return Numeric vector of K distinct positive scales.
#' @export
#' @examples
#' minimax_scales(8)
minimax_scales <- function(K = 22L) {
  K <- as.integer(K)
  if (!(K %in% c(1:11, 14L, 15L, 16L, 17L, 22L))) {
    stop("K must be 1..11 (strict certified), 14..17 (BA-cert), or 22 (certified signed frontier).")
  }
  ## Qualified rather than bare, for the LazyData lookup reason given in
  ## `minimax_error()`.
  LogisNormIntegral::minimax_constants[[paste0("K", K)]]$scales
}

#' Small-v asymptotic series for Pirjol's logistic-normal integral
#'
#' Asymptotic series for the small-variance regime, derived via the
#' transformed Mordell integral representation. Approximates Pirjol's
#' \eqn{\varphi(\eta, v) = 1 - I(\eta, v)}, the reflected form of the integral
#' computed by [logistic_normal_integral].
#'
#' @section Regime of validity:
#' The series is asymptotic, not convergent, so a fixed truncation has a floor
#' that grows with `v`. Against adaptive quadrature at `order = 1` and
#' \eqn{\eta \in [-4, 4]}, the worst-case absolute error is about 1.6e-2 at
#' v = 1, 1.2e-1 at v = 4, and 1.7e-1 at v = 6. Raising `order` helps only while
#' `v` is small. For a bound that holds uniformly in \eqn{(\eta, v)}, use
#' [logistic_normal_integral], whose certified mixture meets
#' \eqn{\Delta_K^\star}. These functions are provided for analytical work in the
#' asymptotic regimes.
#'
#' @param eta Mean of the Gaussian predictor (scalar or vector).
#' @param v Variance of the Gaussian predictor (small, scalar or vector).
#' @param order Number of correction terms (0, 1, or 2)
#' @return Series-truncated approximation to
#'   \eqn{\varphi(\eta, v) = 1 - I(\eta, v)}.
#' @export
#' @examples
#' # Small-variance series, vectorised over eta and v
#' phi_small_v(c(-1, 0, 1), v = 0.5, order = 2)
phi_small_v <- function(eta, v, order = 1) {
  .assert_order(order, "phi_small_v")
  beta <- eta - v / 2
  prefactor <- 0.5 * exp(v / 8 - eta / 2)
  c <- cosh(beta / 2)
  term0 <- 1 / c
  res <- term0
  if (order >= 1) {
    term1 <- 0.25 * (c^2 - 2) / c^3
    res <- res + (v / 2) * term1
  }
  if (order >= 2) {
    # Closed-form 4th derivative of sech(x/2) at x = beta (the heat-kernel
    # term). d^4/dx^4 sech(x/2) = (1/16)(sech - 20 sech^3 + 24 sech^5)(x/2).
    # The large-v branch uses the identical expression in cos. A 5-point
    # finite difference divides cancellation noise by h^4, so the correction
    # was dominated by roundoff (order-2 error worse than order-1).
    f4 <- (1 / 16) * (1 / c - 20 / c^3 + 24 / c^5)
    res <- res + 0.5 * (v / 2)^2 * f4
  }
  prefactor * res
}

# Scalar contract for the large-v and dispatch asymptotic functions. Their
# domain mapping (floor(eta / v), the reflection recursion, the v <= 6
# dispatch) is per-element, unlike phi_small_v which is plain vectorised
# arithmetic. Guard with a clear message instead of a cryptic length error.
.assert_scalar_eta_v <- function(eta, v, fn) {
  if (length(eta) != 1L || length(v) != 1L) {
    stop(fn, " is scalar in `eta` and `v`. Use mapply() or Map() over ",
         "(eta, v) to evaluate on a vector of inputs.")
  }
}

# Highest asymptotic correction the series actually implement. Both phi_small_v
# and phi_large_v_base stop at the `order >= 2` branch, and the .small_terms and
# .large_terms tables carry four entries, exactly enough for .trunc_ratio to read
# the first omitted term at orders 0, 1, and 2. Without this guard order 3 read
# past the end of those tables into an NA comparison, and any higher order
# silently returned the order-2 value.
.MAX_ASYMPTOTIC_ORDER <- 2L

.assert_order <- function(order, fn) {
  if (length(order) != 1L || !is.numeric(order) || !is.finite(order) ||
      order != as.integer(order) || order < 0L ||
      order > .MAX_ASYMPTOTIC_ORDER) {
    stop(fn, " supports `order` 0, 1, or ", .MAX_ASYMPTOTIC_ORDER,
         " only, and got ", paste(format(order), collapse = ", "),
         ". The series implement no higher correction.", call. = FALSE)
  }
}

#' Large-v asymptotic base function on the fundamental domain
#'
#' Evaluates the large-variance series derived from Zwegers' Mordell integral.
#' Valid on the fundamental domain 0 < eta < v. Approximates Pirjol's
#' \eqn{\varphi(\eta, v) = 1 - I(\eta, v)}.
#'
#' @section Regime of validity:
#' The correction terms carry powers of \eqn{2\pi^2 / v}, which is below one
#' only for \eqn{v > 2\pi^2 \approx 19.7}, and the series has poles at the two
#' edges eta = 0 and eta = v of the fundamental domain. It is therefore
#' informative only for large `v` and for eta well inside the domain. Against
#' adaptive quadrature at `order = 1` and \eqn{\eta \in [0.15v, 0.85v]}, the
#' worst-case absolute error is about 6.9 at v = 6, 2.9 at v = 10, 0.83 at
#' v = 20, and 7.7e-2 at v = 60. On \eqn{6 \le v \le 25} neither this series nor
#' [phi_small_v] is accurate. Use [logistic_normal_integral] there, and wherever
#' a bound uniform in \eqn{(\eta, v)} is wanted.
#'
#' @param eta Mean of the Gaussian predictor (scalar).
#' @param v Variance of the Gaussian predictor (scalar).
#' @param order Number of correction terms (0, 1, or 2)
#' @return A numeric scalar, the order-truncated large-variance asymptotic
#'   approximation to \eqn{\varphi(\eta, v) = 1 - I(\eta, v)} on the
#'   fundamental domain 0 < eta < v.
#' @export
#' @examples
#' # Large-variance series on the fundamental domain 0 < eta < v
#' phi_large_v_base(10, v = 20, order = 2)
phi_large_v_base <- function(eta, v, order = 1) {
  .assert_scalar_eta_v(eta, v, "phi_large_v_base")
  .assert_order(order, "phi_large_v_base")
  if (eta <= 0 || eta >= v) stop("phi_large_v_base is only valid for 0 < eta < v")
  a <- pi * (1 - 2 * eta / v)
  prefactor <- sqrt(pi / (2 * v)) * exp(-eta^2 / (2 * v))
  c <- cos(a / 2)
  term0 <- 1 / c
  res <- prefactor * term0
  if (order >= 1) {
    term1 <- - (2 * pi^2 / v) * (0.25 * (2 - c^2) / c^3)
    res <- res + prefactor * term1
  }
  if (order >= 2) {
    term2 <- 0.5 * (2 * pi^2 / v)^2 * (0.125 * (12 / c^5 - 10 / c^3 + 0.5 / c))
    res <- res + prefactor * term2
  }
  res
}

#' Global continuation of the large-v asymptotic series
#'
#' Maps any (eta, v) into the fundamental domain via Pirjol's recurrences. The
#' recurrences are exact, so the accuracy is that of [phi_large_v_base] at the
#' reduced argument and its regime of validity carries over unchanged.
#'
#' @param eta Mean of the Gaussian predictor (scalar).
#' @param v Variance of the Gaussian predictor (scalar).
#' @param order Number of correction terms (0, 1, or 2)
#' @return A numeric scalar, the large-variance asymptotic approximation to
#'   \eqn{\varphi(\eta, v) = 1 - I(\eta, v)} at (eta, v), continued from the
#'   fundamental domain to all real eta by Pirjol's recurrences.
#' @export
#' @examples
#' # Large-variance series with automatic domain continuation
#' phi_large_v_global(10, v = 20, order = 2)
#' phi_large_v_global(-5, v = 20, order = 2)   # reflection for eta < 0
phi_large_v_global <- function(eta, v, order = 1) {
  .assert_scalar_eta_v(eta, v, "phi_large_v_global")
  .assert_order(order, "phi_large_v_global")
  if (eta < 0) return(1 - phi_large_v_global(-eta, v, order))
  n <- floor(eta / v)
  r <- eta - n * v
  if (r == 0) {
    phi_val <- 0.5
  } else {
    phi_val <- phi_large_v_base(r, v, order)
  }
  if (n > 0) {
    for (k in 0:(n - 1)) {
      x <- r + k * v
      phi_val <- exp(-x - v / 2) * (1 - phi_val)
    }
  }
  phi_val
}

## Bracket terms of each series, so a caller can weigh one branch against the
## other. Both are Sum_k T_k times a branch prefactor, and the prefactor cancels
## in |T_{n+1}| / |Sum_{k<=n} T_k|, which makes that ratio comparable across
## branches. Even derivatives of sech and sec, at u = argument/2:
##   sech: S, S - 2S^3, S - 20S^3 + 24S^5, S - 182S^3 + 840S^5 - 720S^7
##   sec : C, 2C^3 - C, 24C^5 - 20C^3 + C, 720C^7 - 840C^5 + 182C^3 - C
.small_terms <- function(eta, v) {
  C <- 1 / cosh((eta - v / 2) / 2)
  g <- c(C, (C - 2 * C^3) / 4, (C - 20 * C^3 + 24 * C^5) / 16,
         (C - 182 * C^3 + 840 * C^5 - 720 * C^7) / 64)
  k <- 0:3
  (v / 2)^k / factorial(k) * g
}
.large_terms <- function(eta, v) {
  C <- 1 / cos(pi * (1 - 2 * eta / v) / 2)
  kap <- c(C, (2 * C^3 - C) / 4, (24 * C^5 - 20 * C^3 + C) / 16,
           (720 * C^7 - 840 * C^5 + 182 * C^3 - C) / 64)
  j <- 0:3
  (-2 * pi^2 / v)^j / factorial(j) * kap
}
## First-omitted-term diagnostic, the standard truncation estimate for an
## asymptotic series. Small means the series is still helping at this order.
.trunc_ratio <- function(terms, order) {
  s <- sum(terms[seq_len(order + 1L)])
  if (!is.finite(s) || s == 0) return(Inf)
  abs(terms[order + 2L] / s)
}
## The large branch is continued from the reduced argument, so judge it there.
.large_ratio <- function(eta, v, order) {
  rr <- abs(eta) - floor(abs(eta) / v) * v
  if (rr == 0) 0 else .trunc_ratio(.large_terms(rr, v), order)
}
## Shared unreliable-truncation warning. `ratios` reports the diagnostics the
## caller was in a position to judge, and `branch` names the series that supplied
## the returned value.
.warn_truncation <- function(eta, v, ratios, branch) {
  warning(sprintf(
    paste0("phi_analytic: the selected asymptotic series is not reliable at (eta = %g, v = %g). ",
           "%s, and the %s-v branch was returned. Use logistic_normal_integral() for a value ",
           "with a certified bound."),
    eta, v, ratios, branch), call. = FALSE)
}

#' Adaptive dispatch between the asymptotic regimes
#'
#' Evaluates both the small-variance and the large-variance series and returns
#' whichever is better behaved at this \eqn{(\eta, v)}, judged by the size of the
#' first omitted term relative to the partial sum. A branch returning a value
#' outside the unit interval is rejected outright when the other does not, since
#' \eqn{\varphi(\eta, v)} is a probability. When both leave the unit interval
#' there is nothing to reject in favour of, so the smaller-ratio branch is
#' returned and the value it returned is warned on.
#'
#' @section Regime of validity:
#' Selection cannot manufacture accuracy that neither series has. The small-v
#' series expands in \eqn{v/2} and the large-v series in \eqn{2\pi^2/v}, so on
#' roughly \eqn{6 \le v \le 25} both are inaccurate and the better of the two may still
#' be wrong in the first or second decimal. Use [logistic_normal_integral] for
#' values. The certified mixture covers that band at \eqn{\Delta_K^\star}, and
#' these series reproduce the paper's analytical forms. Set `warn = FALSE` to
#' silence the warning.
#'
#' @section What the warning does and does not detect:
#' The warning compares the first omitted term of the selected truncation
#' against its own partial sum. That ratio estimates the remainder only while
#' the terms are still decreasing, and inside the band they are not. At
#' \eqn{v \approx 10} the small-variance series expands in \eqn{v/2 \approx 5},
#' so its four tabulated terms grow rather than decay and the ratio ceases to
#' estimate the truncation error, falling below the threshold by virtue of the
#' magnitude of its own denominator.
#' The warning is therefore one-sided. Its presence means the series are
#' visibly failing, and its absence is not evidence of accuracy.
#'
#' Measured on the adaptive selection against adaptive quadrature over
#' \eqn{\eta \in [0.1, 20]} in steps of 0.1 and \eqn{v \in [0.5, 50]} in steps
#' of 0.25 at orders 0, 1, and 2, a grid of 119400 arguments, the warning fires
#' on 62.9 per cent of the grid and the largest error returned without one is
#' 0.30, at \eqn{(\eta, v) = (2.3, 16.5)} and `order = 2`. The range test on the
#' returned value covers the 207 arguments where both branches leave
#' \eqn{[0, 1]}. At 3 of them it is the only signal, a diverged truncation
#' carrying a small ratio whenever its partial sum is large against its first
#' omitted term, and those 3 held the worst error the ratio test alone left
#' unwarned, 0.44 at \eqn{(\eta, v) = (3, 11)} and `order = 2`. No ratio
#' threshold removes what remains. Every ratio criterion tested that leaves no
#' such case unwarned raises a warning on at least 88 per cent of all arguments,
#' including at least 79 per cent of those where the returned value is accurate
#' to better than 0.005, and the error is not monotone in the expansion
#' parameter either. A branch forced through `method` is outside that
#' measurement. It runs the same two tests, and the unit-interval rejection is
#' unavailable, a named branch leaving the selection no alternative to prefer,
#' so the range test on the returned value is the only protection the selection
#' can add. For a value carrying a bound rather than a heuristic,
#' use [logistic_normal_integral].
#'
#' @param eta Mean of the Gaussian predictor (scalar).
#' @param v Variance of the Gaussian predictor (scalar).
#' @param order Number of correction terms (0, 1, or 2)
#' @param method `"adaptive"` (default) selects per argument, `"small"` and
#'   `"large"` force a branch.
#' @param warn Logical. Warn when the first-omitted-term diagnostic of the
#'   selected branch is large at this \eqn{(\eta, v)}, on a branch forced through
#'   `method` as well as on the adaptively selected one, and on either path also
#'   when the value returned falls outside \eqn{[0, 1]}. The diagnostic is a
#'   heuristic rather than a bound, so silence does not certify the value. See
#'   the section below on what it does and does not detect. Default `TRUE`.
#' @return A numeric scalar, the asymptotic-series approximation to
#'   \eqn{\varphi(\eta, v) = 1 - I(\eta, v)}, carrying an attribute `branch`
#'   recording which series supplied it.
#' @export
#' @examples
#' # Selection by argument rather than by a fixed variance threshold
#' phi_analytic(0.5, v = 1, order = 2)
#' phi_analytic(10, v = 20, order = 2)
#' # On the band where neither series works, use the certified mixture instead
#' 1 - logistic_normal_integral(1, 6.5, K = 22)
phi_analytic <- function(eta, v, order = 1,
                         method = c("adaptive", "small", "large"), warn = TRUE) {
  .assert_scalar_eta_v(eta, v, "phi_analytic")
  .assert_order(order, "phi_analytic")
  method <- match.arg(method)

  ## A forced branch runs the same warn logic as the adaptive path. The branch
  ## the caller names is the branch returned, so its own first-omitted-term
  ## ratio is the diagnostic to test, by the selected-branch rule given below.
  ## Returning before this check left the documented `warn` contract inoperative
  ## on both forced branches, which is the silent failure the adaptive path was
  ## corrected for. The ratio is not sufficient on its own here. A diverged
  ## truncation carries a small ratio whenever its partial sum is large against
  ## its first omitted term, so the ratio test is quietest exactly where the
  ## value is worst. At (eta, v) = (3.5, 12), order 2, the forced small branch
  ## returns -0.289 against a true 0.184 at a ratio of 0.056. The adaptive path
  ## catches that case by rejecting a branch outside [0, 1] in favour of the
  ## other one, and a forced branch has no other one, so the range test is
  ## applied to the returned value instead.
  if (method != "adaptive") {
    value  <- if (method == "small") phi_small_v(eta, v, order)
              else phi_large_v_global(eta, v, order)
    r_pick <- if (method == "small") .trunc_ratio(.small_terms(eta, v), order)
              else .large_ratio(eta, v, order)
    out_of_range <- !is.finite(value) || value < 0 || value > 1
    if (isTRUE(warn) && (r_pick > 0.25 || out_of_range)) {
      .warn_truncation(
        eta, v,
        if (out_of_range) {
          sprintf(paste0("The truncation has left the unit interval, returning %g at a ",
                         "first-omitted-term ratio of %.2g"), value, r_pick)
        } else {
          sprintf("The first-omitted-term ratio is %.2g", r_pick)
        },
        method)
    }
    return(structure(value, branch = method))
  }

  small <- phi_small_v(eta, v, order)
  large <- phi_large_v_global(eta, v, order)
  r_small <- .trunc_ratio(.small_terms(eta, v), order)
  r_large <- .large_ratio(eta, v, order)

  ok <- function(x) is.finite(x) && x >= 0 && x <= 1
  pick <- if (ok(small) && !ok(large)) "small"
          else if (ok(large) && !ok(small)) "large"
          else if (r_small <= r_large) "small" else "large"
  value <- if (pick == "small") small else large

  ## Judge the warning on the branch actually returned, not on min(r_small,
  ## r_large). The unit-interval rule above can override the ratio comparison,
  ## so the returned branch is not always the argmin, and the minimum would then
  ## read the rejected branch's ratio. That suppression is anti-correlated with
  ## reliability, since rejection fires when the small branch has diverged out
  ## of [0, 1], and a diverged truncation carries a small ratio whenever its
  ## partial sum is large against its first omitted term. At (eta, v) =
  ## (3.5, 12), order 2, the small branch returned -0.289 and was rejected while
  ## its ratio read 0.056, so the large branch was returned silently at ratio
  ## 4.66 with an error of 0.53 on a quantity confined to [0, 1]. Testing the
  ## selected ratio is strictly conservative, min(...) > 0.25 implying both
  ## ratios exceed 0.25 and therefore implying the selected one does.
  r_pick <- if (pick == "small") r_small else r_large

  ## The unit-interval rule can reject a branch only when the other one is in
  ## range, so it cannot fire when both have diverged. The returned value is
  ## then outside [0, 1] itself, and the ratio test is quietest exactly there,
  ## so the range test is applied to the returned value as well. By the
  ## selection rule above, a returned value outside [0, 1] implies that both
  ## branches were outside it. At (eta, v) = (3, 11), order 2, both leave the
  ## unit interval, the small branch is returned at -0.226 on a ratio of 0.063,
  ## and this test is the only one that fires.
  out_of_range <- !ok(value)

  if (isTRUE(warn) && (r_pick > 0.25 || out_of_range)) {
    .warn_truncation(
      eta, v,
      if (out_of_range) {
        sprintf(paste0("Both branches have left the unit interval, returning %g at ",
                       "first-omitted-term ratios of %.2g (small-v) and %.2g (large-v)"),
                value, r_small, r_large)
      } else {
        sprintf("The first-omitted-term ratios are %.2g (small-v) and %.2g (large-v)",
                r_small, r_large)
      },
      pick)
  }
  structure(value, branch = pick)
}
