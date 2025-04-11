#' Kernel Density Estimation
#'
#' This function estimates the probability density function using Kernel density estimation with Gaussian kernels.
#'
#' @param data A numeric vector of data points.
#' @param points A numeric vector of points at which to estimate the density.
#' @param bandwidth The bandwidth (smoothing parameter) for the Kernel.
#' @return A numeric vector of the estimated density values at the specified points.
#' @examples
#' \dontrun{
#' # Generate sample data
#' data <- rnorm(500, mean = 0, sd = 1)
#' 
#' # Define points at which to estimate density
#' points <- seq(-3, 3, length.out = 100)
#' 
#' # Calculate optimal bandwidth using Silverman's rule of thumb
#' bandwidth <- 1.06 * sd(data) * length(data)^(-1/5)
#' 
#' # Estimate density
#' densities <- kernel_density_estimation(data, points, bandwidth)
#' 
#' # Plot the results
#' plot(points, densities, type = "l", main = "Kernel Density Estimation")
#' }
#' @export
kernel_density_estimation = function(data, points, bandwidth) {
  n_data = length(data)
  n_points = length(points)
  
  kernel_matrix = dnorm((points %*% t(rep(1, n_data)) - t(data %*% t(rep(1, n_points)))) / bandwidth)
  estimated_density = as.vector((kernel_matrix %*% rep(1, n_data)) / (n_data * bandwidth))
  
  return(estimated_density)
}

#' Kernel Distribution Estimation
#'
#' This function estimates the cumulative distribution function (CDF) using Kernel density estimation with Gaussian kernels.
#'
#' @param data A numeric vector of data points.
#' @param points A numeric vector of points at which to estimate the distribution.
#' @param bandwidth The bandwidth (smoothing parameter) for the Kernel.
#' @return A numeric vector of the estimated distribution values at the specified points.
#' @examples
#' \dontrun{
#' # Generate sample data
#' data <- rnorm(500, mean = 0, sd = 1)
#' 
#' # Define points at which to estimate distribution
#' points <- seq(-3, 3, length.out = 100)
#' 
#' # Calculate optimal bandwidth using Silverman's rule of thumb
#' bandwidth <- 1.06 * sd(data) * length(data)^(-1/5)
#' 
#' # Estimate CDF
#' cdf <- kernel_distribution_estimation(data, points, bandwidth)
#' 
#' # Plot the results
#' plot(points, cdf, type = "l", main = "Kernel Distribution Estimation")
#' }
#' @export
kernel_distribution_estimation = function(data, points, bandwidth) {
  n_data = length(data)
  n_points = length(points)
  
  kernel_matrix = pnorm((points %*% t(rep(1, n_data)) - t(data %*% t(rep(1, n_points)))) / bandwidth)
  estimated_distribution = as.vector((kernel_matrix %*% rep(1, n_data)) / n_data)
  
  return(estimated_distribution)
}

#' Evaluate Kernel Estimation at a Point
#'
#' This function evaluates the Kernel estimation at a specific point based on the mesh by interpolating
#' between the closest mesh points.
#'
#' @param point A numeric value for the point at which to evaluate the function.
#' @param function_values A numeric vector of the estimated function values (from Kernel estimation).
#' @param mesh A numeric vector of the mesh/grid used in the estimation.
#' @return The evaluated value at the specified point.
#' @examples
#' \dontrun{
#' # Generate sample data and estimate density
#' data <- rnorm(500, mean = 0, sd = 1)
#' mesh <- seq(-3, 3, length.out = 100)
#' bandwidth <- 1.06 * sd(data) * length(data)^(-1/5)
#' density_values <- kernel_density_estimation(data, mesh, bandwidth)
#' 
#' # Evaluate density at a specific point not in the mesh
#' point <- 1.5
#' evaluated_density <- evaluate_kernel_estimation(point, density_values, mesh)
#' print(evaluated_density)
#' }
#' @export
evaluate_kernel_estimation = function(point, function_values, mesh) {
  sorted_function_values = sort(function_values)
  n_points = length(mesh)
  position = sum((mesh < point))
  r_index = 1
  
  if (position == n_points) {
    r_index = n_points - 1
  }
  
  if (position < n_points && position > 1) {
    r_index = position
  }
  
  value = mean(c(function_values[r_index], function_values[r_index + 1]))
  
  return(value)
}

#' Find Quantile of Kernel Estimation
#'
#' This function finds the quantile (inverse CDF) of the Kernel estimation at a specific probability point
#' by interpolating between the closest function values.
#'
#' @param point A numeric value representing the probability (between 0 and 1) for which to find the quantile.
#' @param function_values A numeric vector of the estimated CDF values (from Kernel estimation).
#' @param mesh A numeric vector of the mesh/grid used in the estimation.
#' @return The quantile value corresponding to the specified probability point.
#' @examples
#' \dontrun{
#' # Generate sample data and estimate CDF
#' data <- rnorm(500, mean = 0, sd = 1)
#' mesh <- seq(-3, 3, length.out = 100)
#' bandwidth <- 1.06 * sd(data) * length(data)^(-1/5)
#' cdf_values <- kernel_distribution_estimation(data, mesh, bandwidth)
#' 
#' # Find the value at the 75th percentile
#' prob <- 0.75
#' quantile_value <- inverse_kernel_estimation(prob, cdf_values, mesh)
#' print(quantile_value)
#' }
#' @export
inverse_kernel_estimation = function(point, function_values, mesh) {
  sorted_function_values = sort(function_values)
  n_points = length(mesh)
  position = sum(function_values < point)
  r_index = 1
  
  if (position == n_points) {
    r_index = n_points - 1
  }
  
  if (position < n_points && position > 1) {
    r_index = position
  }
  
  value = mean(c(mesh[r_index], mesh[r_index + 1]))
  
  return(value)
}

#' Kernel Eta Estimation
#'
#' This function estimates the population eta using Kernel density estimation for controls and cases samples.
#' It calculates the ROC curve and its derivative based on non-parametric kernel density estimates, then
#' uses these to compute the eta value.
#'
#' @param controls A numeric vector of data points representing the controls sample.
#' @param cases A numeric vector of data points representing the cases sample.
#' @param method A character string specifying the bandwidth selection method. Options are 'optimal' (Silverman's rule), 
#'        'iqr' (robust rule based on interquartile range), or 'hscv' (least squares cross-validation from ks package).
#' @param t0 A numeric value indicating the cutoff point for integration. Only mesh points up to t0 are considered. Default is 1.
#' @param mesh_size_kernel The number of mesh points used for the estimation. Default is 1000.
#' @param box_cox A logical value indicating whether to apply a Box-Cox transformation. Default is FALSE.
#' @return A numeric value representing the estimated eta measure.
#' 
#' @importFrom ks hscv
#' 
#' 
#' @examples
#' \dontrun{
#' # Generate sample data with moderate separation
#' controls <- rnorm(1000, mean = 0, sd = 1)
#' cases <- rnorm(1000, mean = 1, sd = 1) 
#' 
#' # Estimate eta using optimal bandwidth selection
#' eta_optimal <- kernel_eta(controls, cases, method = "optimal")
#' print(eta_optimal)
#' 
#' # Estimate eta using IQR-based bandwidth selection
#' eta_iqr <- kernel_eta(controls, cases, method = "iqr")
#' print(eta_iqr)
#' 
#' # Estimate eta with Box-Cox transformation (useful for skewed data)
#' gamma_controls <- rgamma(1000, shape = 2, rate = 1)
#' gamma_cases <- rgamma(1000, shape = 3, rate = 1)
#' eta_boxcox <- kernel_eta(gamma_controls, gamma_cases, method = "optimal", box_cox = TRUE)
#' print(eta_boxcox)
#' 
#' # Using a different cutoff and finer mesh
#' eta_fine <- kernel_eta(controls, cases, method = "optimal", t0 = 0.9, mesh_size_kernel = 2000)
#' print(eta_fine)
#' }
#' 
#' @seealso \code{\link{eta_from_roc_curves}} for the underlying eta calculation,
#'   \code{\link{kernel_density_estimation}} and \code{\link{kernel_distribution_estimation}} for the
#'   kernel estimation functions, \code{\link{apply_box_cox}} for details on the Box-Cox transformation
#' 
#' @export
kernel_eta = function(controls, cases, method, t0 = 1, mesh_size_kernel = 1000, box_cox = FALSE) {
  
  # Apply Box-Cox transformation if specified
  if (box_cox) {
    transformed_data = apply_box_cox(controls, cases)
    controls = transformed_data$transformed_x
    cases = transformed_data$transformed_y
  }
  
  # Combine the samples
  combined_sample = c(controls, cases)
  
  # Bandwidth selection methods
  if (method == 'optimal') {
    bandwidth_controls = 1.06 * sd(controls) * length(controls)^(-1/5)
    bandwidth_cases = 1.06 * sd(cases) * length(cases)^(-1/5)
  }
  else if (method == 'hscv') {
    bandwidth_controls = ks::hscv(controls)
    bandwidth_cases = ks::hscv(cases)
  }
  else if (method == 'iqr') {
    bandwidth_controls = 0.9 * min(sd(controls), (IQR(controls) / 1.34)) * (length(controls)^(-1/5)) 
    bandwidth_cases = 0.9 * min(sd(cases), (IQR(cases) / 1.34)) * (length(cases)^(-1/5))
  }

  
  # Sort the samples
  sorted_controls = sort(controls)
  sorted_cases = sort(cases)
  
  # Create mesh for estimation
  mesh = seq(min(c(controls, cases)),
             max(c(controls, cases)),
             length.out = mesh_size_kernel)
  
  # Estimate distributions and densities using Kernel estimation
  estimated_dist_controls = kernel_distribution_estimation(sorted_controls, mesh, bandwidth_controls)
  estimated_dist_cases = kernel_distribution_estimation(sorted_cases, mesh, bandwidth_cases)
  estimated_density_controls = kernel_density_estimation(sorted_controls, mesh, bandwidth_controls)
  estimated_density_cases = kernel_density_estimation(sorted_cases, mesh, bandwidth_cases)
  
  # Create probability sequence
  p = seq(0.0001, 0.999, length.out = mesh_size_kernel)
  p_opp = 1 - p
  
  # Initialize ROC and ROC' vectors
  roc = numeric(mesh_size_kernel)
  roc_prime = numeric(mesh_size_kernel)
  
  # Compute ROC and ROC' for each point in the mesh
  for (i in 1:mesh_size_kernel) {
    point = p_opp[i]
    
    inv = inverse_kernel_estimation(point, estimated_dist_controls, mesh)
    numerator = evaluate_kernel_estimation(inv, estimated_density_cases, mesh)
    denominator = evaluate_kernel_estimation(inv, estimated_density_controls, mesh)
    
    roc[i] = 1 - evaluate_kernel_estimation(inv, estimated_dist_cases, mesh)
    roc_prime[i] = numerator / denominator
  }
  
  # Return the eta values
  return(eta_from_roc_curves(roc, roc_prime, t0, p))
}
