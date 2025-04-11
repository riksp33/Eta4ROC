
#' Calculate eta value from ROC curves
#'
#' This function calculates the eta value using the given ROC curves and mesh. The eta is computed based on the formula
#' that involves the differences between the ROC curve values and mesh intervals.
#'
#' @param roc A numeric vector of ROC curve values.
#' @param roc_dx A numeric vector of ROC curve derivatives (roc_prima).
#' @param t0 A numeric value indicating the cutoff point for integration. Only mesh points up to t0 are considered.
#' @param mesh A numeric vector representing the mesh/grid used for integration.
#' @return The standardized eta value calculated using the standarization() function.
#' @details The function computes eta by considering two cases: when roc_dx <= 1 and when roc_dx > 1.
#'   For each mesh interval, it calculates a partial eta value and then sums all values up to the t0 threshold.
#'   The final result is standardized using the standarization() function.
#' @examples
#' \dontrun{
#' # Create sample ROC curve values
#' roc <- c(0, 0.2, 0.5, 0.8, 1)
#' roc_dx <- c(0.5, 1.2, 1.5, 1.8, 2.0)
#' mesh <- seq(0, 1, 0.25)
#' 
#' # Calculate eta value with cutoff at 0.75
#' result <- eta_from_roc_curves(roc, roc_dx, 0.75, mesh)
#' print(result)
#' }
#' @export
eta_from_roc_curves = function(roc, roc_dx, t0, mesh) {
  etas = numeric()
  
  # Compute eta for the first element
  etas[1] = ifelse(roc_dx[1] <= 1,
                    ((roc_dx[1] - 1) ^ 2) * mesh[1],
                    ((roc_dx[1] - 1) ^ 2) / roc_dx[1] * roc[1]
  )
  
  # Compute eta for the remaining elements
  for (i in 2:length(mesh)) {
    etas[i] = ifelse(roc_dx[i] <= 1,
                      ((roc_dx[i] - 1) ^ 2) * (mesh[i] - mesh[i-1]),
                      ((roc_dx[i] - 1) ^ 2) / roc_dx[i] * (roc[i] - roc[i-1])
    )
  }
  
  # Sum all eta values and standardize the result
  lim = sum(mesh <= t0)
  eta = sum(etas[1:lim])
  return(standarization(eta))
}


#' Calculate eta analytically for various probability distributions
#'
#' This function computes the eta value using closed-form expressions for different probability distributions
#' (Gaussian, Log-normal, or Gamma) by transforming their parameters into ROC curve values and derivatives.
#'
#' @param param_1_y First parameter of the non-diseased population distribution:
#'   For Gaussian/Log-normal: mean (mu)
#'   For Gamma: shape parameter
#' @param param_2_y Second parameter of the non-diseased population distribution:
#'   For Gaussian/Log-normal: standard deviation (sigma)
#'   For Gamma: rate parameter
#' @param param_1_x First parameter of the diseased population distribution:
#'   For Gaussian/Log-normal: mean (mu)
#'   For Gamma: shape parameter
#' @param param_2_x Second parameter of the diseased population distribution:
#'   For Gaussian/Log-normal: standard deviation (sigma)
#'   For Gamma: rate parameter
#' @param mesh A numeric vector representing the grid points for integration. 
#'   Default is a sequence from 0.00001 to 0.99999 with 10,000 points.
#' @param case A character string specifying the distribution family: "gaussian", "lognormal", or "gamma".
#' @param t0 A numeric value indicating the cutoff point for integration. Default is 1.
#' @return The calculated eta value after standardization.
#' @details The function works by:
#'   1. Setting up appropriate quantile, density, and cumulative distribution functions for the selected distribution
#'   2. Computing inverse values, ROC curve points, and derivatives based on the distribution parameters
#'   3. Calling eta_from_roc_curves() to calculate the final eta value
#' @examples
#' \dontrun{
#' # Gaussian case with different means, same standard deviation
#' normal_eta <- analytical_eta(0, 1, 2, 1, case = "gaussian")
#' print(normal_eta)
#' 
#' # Log-normal case
#' lognormal_eta <- analytical_eta(0, 0.5, 0.5, 0.75, case = "lognormal")
#' print(lognormal_eta)
#' 
#' # Gamma case with different shape and rate parameters
#' gamma_eta <- analytical_eta(2, 1, 3, 2, case = "gamma", t0 = 0.9)
#' print(gamma_eta)
#' 
#' # Using a custom mesh with fewer points
#' custom_mesh <- seq(0.001, 0.999, length.out = 1000)
#' eta_custom <- analytical_eta(0, 1, 1.5, 1.2, mesh = custom_mesh, case = "gaussian")
#' print(eta_custom)
#' }
#' @seealso \code{\link{eta_from_roc_curves}} for the underlying computation method
#' @export
analytical_eta = function(param_1_y, param_2_y, param_1_x, param_2_x, mesh = seq(0.00001, 0.99999, length.out = 10000), case = "gaussian", t0 = 1) {


  if(case == "gaussian") {
    inverse_dist = function(p, mu, sigma) {qnorm(p, mu, sigma)}
    density_func = function(p, mu, sigma) {dnorm(p, mu, sigma)}
    distribution = function(p, mu, sigma) {pnorm(p, mu, sigma)}
  } else if(case == "lognormal") {
    inverse_dist = function(p, mu, sigma) {qlnorm(p, mu, sigma)}
    density_func = function(p, mu, sigma) {dlnorm(p, mu, sigma)}
    distribution = function(p, mu, sigma) {plnorm(p, mu, sigma)}
  } else if(case == "gamma") {
    inverse_dist = function(p, shape, rate) {qgamma(p, shape = shape, rate = rate)}
    density_func = function(p, shape, rate) {dgamma(p, shape = shape, rate = rate)}
    distribution = function(p, shape, rate) {pgamma(p, shape = shape, rate = rate)}
  } else {
    stop("Non valid input, the param case must be gaussian, lognormal or gamma")
  }
  

    inv = inverse_dist(1-mesh, param_1_y, param_2_y)

    num = density_func(inv, param_1_x, param_2_x)
    denom = density_func(inv, param_1_y, param_2_y)

    roc = 1 - distribution(inv, param_1_x, param_2_x)
    roc_prima = num / denom

    return(eta_from_roc_curves(roc, roc_prima, t0, mesh))
}