#' Half-normal distribution around known survival values.
#'
#' @param survival_mean Numeric value. Mean desired survival.
#' @param survival_min Numeric value. Minimum survival.
#' @param survival_max Numeric value. Maximum survival.
#' @param n_samples Number of samples. Default 1000.
#' @param sigma_start Positive valued scale parameter.
#' @param optimize_sigma Logical. Optimize sigma value selection. Default TRUE.
#'
#' @return Creates a vector with survival values.
#' @importFrom extraDistr qhnorm rhnorm
#' @importFrom scales rescale
#' @export
#'
#' @examples
#' \dontrun{
#' survival_values_hn <- make_half_normal(survival_mean = 0.95244,
#' survival_min = 0.9215,
#' survival_max = 0.98153)
#' hist(survival_values)
#' mean(survival_values)
#' }
make_half_normal <- function(survival_mean = NA, survival_min = NA,
                             survival_max = NA, n_samples = 1000,
                             sigma_start = NA, optimize_sigma = TRUE) {
  # Input validation
  if (is.na(survival_mean) || is.na(survival_min) || is.na(survival_max)) {
    stop("Please provide values for survival_mean, survival_min, and survival_max.")
  }
  if (survival_min >= survival_max || survival_mean < survival_min || survival_mean > survival_max) {
    stop("Invalid input: survival_min must be less than survival_max, and survival_mean must be within this range.")
  }

  mean_survival_ajusted <- survival_mean + (survival_mean * 0.05)
  if (mean_survival_ajusted > 0.99){
    mean_survival_ajusted <- 0.99
  }
  # Estimate sigma parameter of the half-normal distribution.
  # This estimated sigma ensures that when we generate random samples
  # from the half-normal distribution with this sigma, the mean of
  # those samples will be close to the desired mean_survival.
  if (optimize_sigma) {
    # Objective function to minimize the difference between desired and actual mean
    objective_function <- function(sigma) {
      hn <- extraDistr::rhnorm(n_samples, sigma)
      # rescale half-normal to desired range
      sv <- scales::rescale(hn, to = c(survival_min, survival_max))
      # minimize difference
      abs(mean(sv, na.rm = TRUE) - mean_survival_ajusted)
    }

    # Optimize to find the best sigma
    if (is.na(sigma_start)) {
      sigma_start <- sqrt(pi/2) / extraDistr::qhnorm(mean_survival_ajusted, lower.tail = FALSE) # Initial guess if not provided
    }
    # Set reasonable finite bounds for sigma based on the desired mean and range
    lower_bound_sigma <- 0.01 # Small positive value to avoid zero
    upper_bound_sigma <- 10 * sigma_start # Allow for a wide range of exploration
    optimization_result <- optim(sigma_start, objective_function,
                                 method = "Brent",
                                 lower = lower_bound_sigma, upper = upper_bound_sigma)
    sigma <- optimization_result$par
  } else {
    # Without optimizing, use the provided sigma or the default estimation if not provided
    sigma <- ifelse(is.na(sigma_start), sqrt(pi/2) / extraDistr::qhnorm(mean_survival_ajusted, lower.tail = FALSE), sigma_start)
  }

  # Generate half-normal samples
  hn_values <- extraDistr::rhnorm(n_samples, sigma)

  # Rescale survival values within the desired range
  survival_values <- scales::rescale(hn_values, to = c(survival_min, survival_max))

  # Create output dataframe with original values
  # Make sure vectors have equal length.
  sq_range <- seq(max(length(hn_values), length(survival_values)))
  dfout <- data.frame(hn_values = hn_values,
                      survival_values = survival_values[sq_range]
  )
  return(dfout)
}
