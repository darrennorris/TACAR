#' Beta distribution around known survival values.
#'
#' @param survival_mean Numeric value. Mean desired survival.
#' @param survival_min Numeric value. Minimum survival.
#' @param survival_max Numeric value. Maximum survival.
#' @param alpha_start Numeric value. Alpha parameter, the first shape parameter of the beta distribution.
#' @param beta_start Numeric value. Beta parameter, the second shape parameter of the beta distribution.
#'
#' @return Creates a data.frame with survival values and beta distribution parameters".
#' @importFrom stats optim rbeta
#' @importFrom scales rescale
#' @export
#'
#' @examples
#' \dontrun{
#' # Here symmetric mean = 30 Y, min = 20 Y and max = 70 Y.
#' adult_survival_symmetric <- TACAR::make_beta(survival_mean = 0.95244,
#'                                      survival_min = 0.9215,
#'                                       survival_max = 0.98153,
#'                                       alpha_start = 200,
#'                                      beta_start = 200)
#' }
make_beta <- function(survival_mean = NA, survival_min = NA,
                      survival_max = NA,
                      alpha_start = 40, beta_start = 2
){

  # Input validation
  if (is.na(survival_mean) || is.na(survival_min) || is.na(survival_max)) {
    stop("Please provide values for survival_mean, survival_min, and survival_max.")
  }
  if (survival_min >= survival_max || survival_mean < survival_min || survival_mean > survival_max) {
    stop("Invalid input: survival_min must be less than survival_max, and survival_mean must be within this range.")
  }


  # Set initial parameter values and constraints.
  # Initial guess for alpha and beta
  initial_params <- c(alpha_start, beta_start)
  # Ensure alpha and beta are positive
  lower_bounds <- c(1, 1)
  upper_bounds <- c(Inf, Inf)

  # Minimize the difference between desired and actual mean
  objective_function_normal <- function(params) {
    p_alpha <- params[1]
    p_beta <- params[2]
    mean_beta <- p_alpha / (p_alpha + p_beta)
    # Minimize the difference between desired and actual mean
    abs(mean_beta - survival_mean)
  }

  objective_function_alpha_less_than_beta <- function(params) {
    p_alpha <- params[1]
    p_beta <- params[2]
    mean_beta <- p_alpha / (p_alpha + p_beta)
    penalty <- ifelse(p_alpha >= p_beta, 1000, 0)
    abs(mean_beta - survival_mean) + penalty
  }

  # Choose the appropriate objective function based on alpha_start and beta_start
  if (alpha_start < beta_start) {
    objective_function <- objective_function_alpha_less_than_beta
  } else {
    objective_function <- objective_function_normal
  }

  # Optimize to find suitable alpha and beta
  optimization_result <- optim(initial_params, objective_function,
                               lower = lower_bounds, upper = upper_bounds,
                               method = "L-BFGS-B")
  # Extract optimized alpha and beta
  alpha_optimized <- optimization_result$par[1]
  beta_optimized <- optimization_result$par[2]

  # Generate random survival values using the optimized parameters
  sv <- rbeta(1000, alpha_optimized, beta_optimized)
  # Keep survival values within the desired range
  all_within_range <- sv[sv >= survival_min & sv <= survival_max]
  # Return data.frame.
  prop_within_range <- length(all_within_range) / length(sv)
  # Check if less than 50% are within range
  if (prop_within_range < 0.5) {
    warning("Less than 50% of survival values are within the desired range. Rescaling values.")

    # Rescale survival values to fit within the range
    within_range_rescaled <- scales::rescale(sv, to = c(survival_min, survival_max))

    # Create output data.frame with rescaled values
    # Make sure vectors have equal length.
    sq_range <- seq(max(length(sv), length(within_range_rescaled)))
    sq_opti <- seq(max(length(sv), length(alpha_optimized)))
    dfout <- data.frame(
      survival_values = sv,
      within_range = within_range_rescaled[sq_range],
      alpha = alpha_optimized,
      beta = beta_optimized
    )
  } else {
    # Create output dataframe with original values
    # Make sure vectors have equal length.
    sq_range <- seq(max(length(sv), length(all_within_range)))
    sq_opti <- seq(max(length(sv), length(alpha_optimized)))
    dfout <- data.frame(survival_values = sv,
                        within_range = all_within_range[sq_range],
                        alpha = alpha_optimized[sq_opti],
                        beta = beta_optimized[sq_opti])
  }

  return(dfout)

}
