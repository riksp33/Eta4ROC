#' Kernel Density Estimation
#'
#' This function estimates the probability density function using Kernel density estimation.
#'
#' @param data A numeric vector of data points.
#' @param points A numeric vector of points at which to estimate the density.
#' @param bandwidth The bandwidth (smoothing parameter) for the Kernel.
#' @return A numeric vector of the estimated density values at the specified points.
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
#' This function estimates the cumulative distribution function (CDF) using Kernel density estimation.
#'
#' @param data A numeric vector of data points.
#' @param points A numeric vector of points at which to estimate the distribution.
#' @param bandwidth The bandwidth (smoothing parameter) for the Kernel.
#' @return A numeric vector of the estimated distribution values at the specified points.
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
#' This function evaluates the Kernel estimation at a specific point based on the mesh.
#'
#' @param point A numeric value for the point at which to evaluate the function.
#' @param function_values A numeric vector of the estimated function values (from Kernel estimation).
#' @param mesh A numeric vector of the mesh/grid used in the estimation.
#' @return The evaluated value at the specified point.
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
#' This function finds the quantile of the Kernel estimation at a specific point.
#'
#' @param point A numeric value representing the quantile to find.
#' @param function_values A numeric vector of the estimated function values (from Kernel estimation).
#' @param mesh A numeric vector of the mesh/grid used in the estimation.
#' @return The quantile value corresponding to the specified point in the mesh.
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
#' This function estimates the population eta using Kernel density estimation for healthy and sick samples.
#' It can apply a Box-Cox transformation if specified.
#'
#' @param controls A numeric vector of data points representing the healthy sample.
#' @param cases A numeric vector of data points representing the sick sample.
#' @param method A character string specifying the bandwidth selection method. Options are 'optimal' or 'hscv'. The latter coming from the ks package.
#' @param mesh_size_kernel The number of mesh points used for the estimation. Default is 1000.
#' @param box_cox A logical value indicating whether to apply a Box-Cox transformation. Default is FALSE.
#' @return A numeric vector of the estimated eta values at the specified points in the mesh.
#' 
#' @importFrom ks hscv
#' 
#' @examples
#' # Example usage without Box-Cox transformation
#' controls = rnorm(1000, mean = 0, sd = 1)
#' cases = rnorm(1000, mean = 1, sd = 1)
#' eta = kernel_eta(controls, cases)
#'
#' # Example usage with Box-Cox transformation
#' eta_boxcox = kernel_eta(controls, cases, box_cox = TRUE)
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
    bandwidth = 1.06 * sd(combined_sample) * length(combined_sample)^(-1/5)
    bandwidth_healthy = 1.06 * sd(controls) * length(controls)^(-1/5)
    bandwidth_sick = 1.06 * sd(cases) * length(cases)^(-1/5)
  }
  if (method == 'hscv') {
    bandwidth = ks::hscv(combined_sample)
    bandwidth_healthy = ks::hscv(controls)
    bandwidth_sick = ks::hscv(cases)
  }
  
  # Sort the samples
  sorted_healthy = sort(controls)
  sorted_sick = sort(cases)
  
  # Create mesh for estimation
  mesh = seq(min(c(controls, cases)),
             max(c(controls, cases)),
             length.out = mesh_size_kernel)
  
  # Estimate distributions and densities using Kernel estimation
  estimated_dist_healthy = kernel_distribution_estimation(sorted_healthy, mesh, bandwidth)
  estimated_dist_sick = kernel_distribution_estimation(sorted_sick, mesh, bandwidth)
  estimated_density_healthy = kernel_density_estimation(sorted_healthy, mesh, bandwidth)
  estimated_density_sick = kernel_density_estimation(sorted_sick, mesh, bandwidth)
  
  # Create probability sequence
  p = seq(0.0001, 0.999, length.out = mesh_size_kernel)
  p_opp = 1 - p
  
  # Initialize ROC and ROC' vectors
  roc = numeric(mesh_size_kernel)
  roc_prime = numeric(mesh_size_kernel)
  
  # Compute ROC and ROC' for each point in the mesh
  for (i in 1:mesh_size_kernel) {
    point = p_opp[i]
    
    inv = inverse_kernel_estimation(point, estimated_dist_healthy, mesh)
    numerator = evaluate_kernel_estimation(inv, estimated_density_sick, mesh)
    denominator = evaluate_kernel_estimation(inv, estimated_density_healthy, mesh)
    
    roc[i] = 1 - evaluate_kernel_estimation(inv, estimated_dist_sick, mesh)
    roc_prime[i] = numerator / denominator
  }
  
  # Return the eta values
  return(eta_from_roc_curves(roc, roc_prime, t0, p))
}
