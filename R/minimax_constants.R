#' Certified minimax constants for the logistic-normal integral
#'
#' @description
#' A list of length 16. Elements `K1, ..., K11` (strict positive certs) each have three fields:
#' \itemize{
#'   \item `K`: integer, the mixture order.
#'   \item `weights`: numeric vector of length K (sums to 1, all positive).
#'   \item `scales`: numeric vector of length K (distinct, positive).
#' }
#'
#' For K = 1, ..., 11, these are the certified global minimax of the
#' approximation problem \eqn{\inf_{f \in F_K} \sup_z |\Lambda(z) - f(z)|}, proved to
#' equioscillate and to be the unique global minimiser (via varisolvency), with the
#' results set out in the package vignette.
#' Source: 100-digit Newton equioscillation
#' solves (`mp.dps = 115`) for K = 1..10 and Sturm-seeded CMA-ES + chain LM +
#' sequential 200k LM in 100-digit arithmetic (Python mpmath) for K = 11. The
#' K = 1..8 solves reproduce the published values of Monahan & Stefanski
#' (Balakrishnan 1992 Table 18.4.1) to about eight significant figures. The
#' K = 9 and K = 10 solves are seeded by CMA-ES + Remez and by
#' perturbed-reference iterated Remez respectively.
#'
#' Element `K22` is the certified \emph{signed} frontier (best certified bound,
#' \eqn{\Delta \le 3.605 \times 10^{-15}}). It adds fields `signed = TRUE`,
#' `sum_abs_weights` (\eqn{\sum_i |w_i| \approx 1.024}, the evaluation-conditioning number),
#' and `cert_form`, and its `weights` include one negative atom while \eqn{\sum_i w_i = 1} still
#' holds.
#'
#' Elements `K14` to `K17` are the relaxed Bernstein--Akhiezer (BA) certs of the positive class,
#' all weights > 0, with K14, K15, and K16 BA-(K+1) and K17 BA-K. They add `signed = FALSE` and
#' `cert_form`.
#'
#' @format List of 16: `K1`--`K11` (strict positive certs), `K14`--`K17` (relaxed BA-certs), and
#'   `K22` (signed frontier).
#' @seealso [minimax_delta], [phi_minimax_mixture], [logistic_normal_integral]
"minimax_constants"

#' Certified minimax error bounds Delta_K*
#'
#' @description
#' Named numeric vector (`K1, ..., K11, K14, ..., K17, K22`). Element `K` is
#' \eqn{\Delta_K^\star}, the
#' certified error bound for the K-term mixture (the supremum-norm error the mixture
#' is verified to meet). See `inst/certs/README.md` for the certification standard and the
#' package vignette for the equioscillation, local-minimality, and global-optimality results.
#'
#' @format Named numeric vector of length 16 (`K1, ..., K11, K14, ..., K17, K22`).
#' @seealso [minimax_constants], [minimax_error]
"minimax_delta"
