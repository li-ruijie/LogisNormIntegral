#' Predict Causal Mediation Effects via Fast Quasi-Bayesian Monte Carlo
#'
#' @description
#' Computes the Natural Indirect Effect (NIE) and Expected Potential Outcomes
#' \eqn{E[Y(r_Y, M(r_M))]} for a binary outcome and continuous mediator.
#'
#' **Rationale for Quasi-Bayesian over the Delta Method:**
#' While the Delta Method can provide analytic confidence intervals for single
#' potential outcomes, it forces a symmetric Normal approximation
#' onto the sampling variance. Indirect effects (like the NIE) are the product
#' of two random variables (the \eqn{\alpha\beta} path). Their true sampling distribution
#' is highly skewed and asymmetric. Applying the Delta Method (Sobel Test) here
#' produces incorrect, artificially symmetric confidence intervals that fail to
#' capture the true uncertainty skew.
#'
#' Quasi-Bayesian Monte Carlo simulation avoids that approximation, at the cost of
#' evaluating the logistic-normal integral by numerical quadrature inside the draw
#' loop, which has historically made it slow. Replacing that quadrature with the
#' certified minimax normal scale mixture (\code{\link{phi_minimax_mixture}},
#' evaluated at the signed frontier K = 22, the lowest certified error
#' 3.60e-15) reduces the inner integral to a closed-form sum accurate to
#' \eqn{\Delta_K^\star}, so tens of thousands of draws complete in well under a
#' second and the skew and asymmetry of the indirect effect are retained in full.
#'
#' @param outcome_model A fitted binomial GLM predicting Y.
#' @param mediator_model A fitted linear model predicting the continuous mediator M.
#' @param treat_val Numeric value of the treatment condition (e.g., 1).
#' @param control_val Numeric value of the control condition (e.g., 0).
#' @param covar_vals A named list or single-row data.frame giving one covariate
#'   profile to condition on. Multi-row input is rejected. To aggregate across
#'   individuals, call once per profile and average the draws (see the group
#'   g-computation example).
#' @param treatment_var Character string naming the treatment variable.
#' @param mediator_var Character string naming the mediator variable.
#' @param n_sim Number of Quasi-Bayesian Monte Carlo simulations (default 10000).
#' @param alpha Significance level for confidence intervals (default 0.05).
#' @param return_draws Logical. If \code{TRUE}, returns the raw Quasi-Bayesian
#'   Monte Carlo simulation draws (default \code{FALSE}).
#'
#' @return A list containing the estimates and asymmetric confidence intervals
#'         for the Natural Indirect Effect (NIE) and the potential outcomes.
#'         If \code{return_draws} is \code{TRUE}, a \code{draws} list is also included.
#'
#' @examples
#' \donttest{
#' # 1. Simulate data
#' set.seed(42)
#' N <- 500
#' treat <- rbinom(N, 1, 0.5)
#' covar <- rnorm(N)
#' med <- 0.5 + 1.2 * treat + 0.3 * covar + rnorm(N, sd = 0.8)
#' logit_p <- -1.0 + 0.5 * treat + 1.5 * med - 0.2 * covar
#' y <- rbinom(N, 1, plogis(logit_p))
#' df <- data.frame(y = y, treat = treat, med = med, covar = covar)
#'
#' # 2. Fit models
#' mod_m <- lm(med ~ treat + covar, data = df)
#' mod_y <- glm(y ~ treat + med + covar, data = df, family = binomial)
#'
#' # ---------------------------------------------------------
#' # Example A: Individual CI (Predicted Probabilities)
#' # ---------------------------------------------------------
#' # Compute for a specific individual profile (covar = 0)
#' res_indiv <- predict_indirect_effect_qb(
#'   outcome_model = mod_y,
#'   mediator_model = mod_m,
#'   covar_vals = data.frame(covar = 0),
#'   treatment_var = "treat",
#'   mediator_var = "med",
#'   return_draws = TRUE
#' )
#' print(res_indiv$NIE)
#'
#' # ---------------------------------------------------------
#' # Example B: Individual CI (Odds Ratio)
#' # ---------------------------------------------------------
#' # Using the same individual simulation draws, calculate the Odds Ratio (OR) distribution
#' odds_11 <- res_indiv$draws$Y_1_M_1 / (1 - res_indiv$draws$Y_1_M_1)
#' odds_10 <- res_indiv$draws$Y_1_M_0 / (1 - res_indiv$draws$Y_1_M_0)
#' OR_NIE <- odds_11 / odds_10
#'
#' cat("Individual OR Estimate:", mean(OR_NIE), "\n")
#' cat("Individual OR 95% CI:", quantile(OR_NIE, c(0.025, 0.975)), "\n")
#'
#' # ---------------------------------------------------------
#' # Example C: Group CI (G-computation across the sample)
#' # ---------------------------------------------------------
#' # To compute the group expected potential outcomes without bias, evaluate the
#' # integral for every individual and average them per simulation draw.
#' n_sim <- 10000
#' group_Y_1_M_1 <- numeric(n_sim)
#' group_Y_1_M_0 <- numeric(n_sim)
#'
#' # Loop over a subset of individuals (first 50) for illustration
#' for (i in 1:50) {
#'   res_i <- predict_indirect_effect_qb(
#'     outcome_model = mod_y,
#'     mediator_model = mod_m,
#'     covar_vals = df[i, , drop = FALSE],
#'     treatment_var = "treat",
#'     mediator_var = "med",
#'     return_draws = TRUE
#'   )
#'   # Accumulate draws
#'   group_Y_1_M_1 <- group_Y_1_M_1 + res_i$draws$Y_1_M_1
#'   group_Y_1_M_0 <- group_Y_1_M_0 + res_i$draws$Y_1_M_0
#' }
#'
#' # Average across the individuals
#' group_Y_1_M_1 <- group_Y_1_M_1 / 50
#' group_Y_1_M_0 <- group_Y_1_M_0 / 50
#' group_NIE <- group_Y_1_M_1 - group_Y_1_M_0
#'
#' cat("Group NIE Estimate:", mean(group_NIE), "\n")
#' cat("Group NIE 95% CI:", quantile(group_NIE, c(0.025, 0.975)), "\n")
#' }
#' @importFrom stats coef vcov quantile delete.response terms model.matrix rchisq
#' @export
predict_indirect_effect_qb <- function(outcome_model, mediator_model,
                                       treat_val = 1, control_val = 0,
                                       covar_vals = list(),
                                       treatment_var = "treat",
                                       mediator_var = "med",
                                       n_sim = 10000,
                                       alpha = 0.05,
                                       return_draws = FALSE) {

  # 1. Quasi-Bayesian Parameter Draws
  #
  # Draws use a Cholesky factor rather than an eigendecomposition. The
  # Cholesky factor of a positive definite matrix is unique (upper triangular
  # with positive diagonal), so a given seed yields the same draws on every
  # BLAS and LAPACK build. MASS::mvrnorm, used here in earlier drafts,
  # factors with eigen(), whose eigenvector signs are not fixed by the matrix,
  # so the same seed produced different (equally valid) draws under a different
  # BLAS and the returned interval moved in its third decimal. The pivoted
  # branch keeps the tolerance of a numerically singular covariance that
  # mvrnorm obtained from its pmax(ev, 0), and is also sign-determinate.
  rmvnorm_chol <- function(n, mu, Sigma) {
    p <- length(mu)
    R <- tryCatch(chol(Sigma), error = function(e) NULL)
    if (is.null(R)) {
      pc <- chol(Sigma, pivot = TRUE)
      R  <- pc[, order(attr(pc, "pivot")), drop = FALSE]
    }
    out <- matrix(stats::rnorm(n * p), nrow = n, ncol = p) %*% R
    out <- sweep(out, 2L, mu, "+")
    dimnames(out) <- list(NULL, names(mu))
    out
  }

  # This routine conditions on a single covariate profile. Multi-row input
  # would be silently reduced to its first row by the model.matrix [1, ]
  # extractions below, so reject it explicitly. To aggregate across
  # individuals, call once per profile and average the draws (see the group
  # g-computation example).
  n_profile <- nrow(as.data.frame(covar_vals))
  if (n_profile > 1L) {
    stop("`covar_vals` must specify a single covariate profile (one row), got ",
         n_profile, " rows. Call once per profile and average the draws to ",
         "aggregate across individuals.")
  }

  sim_theta <- rmvnorm_chol(n_sim, coef(outcome_model), vcov(outcome_model))
  sim_gamma <- rmvnorm_chol(n_sim, coef(mediator_model), vcov(mediator_model))

  # Simulate the residual variance of the mediator
  s_sq_est <- summary(mediator_model)[["sigma"]]^2
  df_resid <- mediator_model[["df.residual"]]
  sim_s_sq <- df_resid * s_sq_est / rchisq(n_sim, df_resid)

  # 2. Build Model Matrices for the requested covariate profile
  df_base <- as.data.frame(covar_vals)
  if (nrow(df_base) == 0) df_base <- data.frame(placeholder = 1)

  # Predictor matrices for M(treat) and M(control)
  df_M1 <- df_base; df_M1[[treatment_var]] <- treat_val
  df_M0 <- df_base; df_M0[[treatment_var]] <- control_val

  X_M1 <- model.matrix(delete.response(terms(mediator_model)), data = df_M1)[1, ]
  X_M0 <- model.matrix(delete.response(terms(mediator_model)), data = df_M0)[1, ]

  # Predictor matrix for Y(treat) excluding the mediator. The natural indirect
  # effect holds the exposure at treat_val in the outcome model for both terms,
  # so no control-exposure outcome row is needed.
  df_Y1 <- df_base; df_Y1[[treatment_var]] <- treat_val; df_Y1[[mediator_var]] <- 0

  X_Y1 <- model.matrix(delete.response(terms(outcome_model)), data = df_Y1)[1, ]

  # Matrix to extract the 'b' coefficient (the mediator effect on Y)
  df_Y_med1 <- df_Y1; df_Y_med1[[mediator_var]] <- 1
  X_Y_med1 <- model.matrix(delete.response(terms(outcome_model)), data = df_Y_med1)[1, ]
  b_vec <- X_Y_med1 - X_Y1

  # 3. Vectorised certified minimax scale mixture across all n_sim draws.
  # Delegates to the package's certified K = 22 signed-frontier mixture (the lowest
  # certified error), so the deployed inner integral inherits the uniform bound
  # |I(eta, v) - I_hat| <= Delta_22, about 3.6e-15 (see logistic_normal_integral,
  # minimax_error(22)). I(eta, v) = E[Lambda(Y)] = sum_i w_i Phi(s_i eta / sqrt(1 + s_i^2 v)).
  minimax_phi_vec <- function(eta, v) {
    logistic_normal_integral(eta, v, K = 22L)
  }

  # 4. Compute Linear Predictors across all draws
  mu_1 <- as.numeric(sim_gamma %*% X_M1)
  mu_0 <- as.numeric(sim_gamma %*% X_M0)

  a_1 <- as.numeric(sim_theta %*% X_Y1)
  b   <- as.numeric(sim_theta %*% b_vec)

  # 5. Compute Eta and V for the potential outcomes under Treatment
  # E[Y(1, M(1))]
  eta_11 <- a_1 + b * mu_1
  v_11   <- (b^2) * sim_s_sq

  # E[Y(1, M(0))]
  eta_10 <- a_1 + b * mu_0
  v_10   <- (b^2) * sim_s_sq

  # 6. Evaluate the integral for all n_sim draws in one vectorised call
  Y_1_M_1 <- minimax_phi_vec(eta_11, v_11)
  Y_1_M_0 <- minimax_phi_vec(eta_10, v_10)

  # 7. Compute the Natural Indirect Effect
  NIE <- Y_1_M_1 - Y_1_M_0

  # Extract properly asymmetric confidence intervals
  res <- list(
    NIE = list(
      estimate = mean(NIE),
      ci_lower = unname(quantile(NIE, alpha / 2)),
      ci_upper = unname(quantile(NIE, 1 - alpha / 2))
    ),
    Y_1_M_1 = list(
      estimate = mean(Y_1_M_1),
      ci_lower = unname(quantile(Y_1_M_1, alpha / 2)),
      ci_upper = unname(quantile(Y_1_M_1, 1 - alpha / 2))
    ),
    Y_1_M_0 = list(
      estimate = mean(Y_1_M_0),
      ci_lower = unname(quantile(Y_1_M_0, alpha / 2)),
      ci_upper = unname(quantile(Y_1_M_0, 1 - alpha / 2))
    )
  )

  if (return_draws) {
    res <- `[[<-`(res, "draws", list(
      NIE = NIE,
      Y_1_M_1 = Y_1_M_1,
      Y_1_M_0 = Y_1_M_0
    ))
  }

  return(res)
}
