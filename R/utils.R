# General utils used trough the main functions

#' Calculate the bias of the estimator
#'
#' @param estimates A numeric vector of estimated values.
#' @param true_value The true value of the estimator.
#' @return The bias of the estimator (difference between mean of estimates and true value).
#' @export
get_bias = function(estimates, true_value) {
  bias = mean(estimates) - true_value
  return(bias)
}


#' Calculate the Root Mean Square Error (RMSE)
#'
#' @param estimates A numeric vector of estimated values.
#' @param true_value The true value of the estimator.
#' @return The RMSE, a measure of the difference between estimated and true values.
#' @export
get_rmse = function(estimates, true_value) {
  rmse = sqrt(mean((estimates - true_value)^2))
  return(rmse)
}

#' Calculate mux value based on (6) from Faraggi, D., & Reiser, B. (2002). Estimation of the area under the ROC curve. Statistics in medicine, 21(20), 3093-3106.
#'
#' @param auc The area under the curve (AUC) value.
#' @param sigma_x Standard deviation of X.
#' @param sigma_y Standard deviation of Y.
#' @param mu_y The mean of Y.
#' @return The calculated mux value.
#' @export
get_mux = function(auc, sigma_x, sigma_y, mu_y) {
  mux = sqrt(sigma_y^2 + sigma_x^2) * qnorm(auc) + mu_y
  return(mux)
}

standarization = function(x) {
  return(log(x +1)/ (1 + log(x+1)) )
}


#' Perform Box-Cox transformation on two variables
#'
#' This function applies the Box-Cox transformation to two numeric vectors. If any values are less than or equal to zero,
#' a small constant is added to both vectors to make them positive.
#'
#' @param x A numeric vector (first variable).
#' @param y A numeric vector (second variable).
#' @param print_lambda A logical value; if TRUE, the lambda value (Box-Cox transformation parameter) is printed.
#' @return A list containing the transformed vectors [transformed_x, transformed_y].
#' @export
apply_box_cox = function(x, y, print_lambda = FALSE) {
  # Count negative observations
  num_negative_obs = sum(x <= 0) + sum(y <= 0)
  
  # If there are negative or zero values, add a small constant to both x and y
  if (num_negative_obs > 0) {
    constant = -min(min(x), min(y)) + 5e-4
    x = x + constant
    y = y + constant
  }
  
  # Likelihood function for Box-Cox
  likbox = function(h, data, n) {
    m = length(data) - n
    x_data = data[1:n]
    y_data = data[(n + 1):length(data)]
    
    # Calculate transformed data based on h
    if (abs(h) < 1e-5) {
      xh = log(x_data)
      yh = log(y_data)
    } else {
      xh = ((x_data^h) - 1) / h
      yh = ((y_data^h) - 1) / h
    }
    
    # Log-likelihood function for Box-Cox
    oout = -n / 2 * log(sum((xh - sum(xh) / n)^2) / n) -
      m / 2 * log(sum((yh - sum(yh) / m)^2) / m) + 
      (h - 1) * (sum(log(x_data)) + sum(log(y_data)))
    
    return(-oout)
  }
  
  # Initial guess for lambda
  h_ini = -0.6
  hhat = optim(h_ini, likbox, data = c(x, y), n = length(x), method = "BFGS")$par
  
  # Apply Box-Cox transformation
  if (abs(hhat) < 1e-5) {
    transformed_x = log(x)
    transformed_y = log(y)
  } else {
    transformed_x = ((x^hhat) - 1) / hhat
    transformed_y = ((y^hhat) - 1) / hhat
  }
  
  # Print the lambda value if requested
  if (print_lambda) {
    cat('The lambda value for the transformation is: ', hhat, '\n')
  }
  
  # Return the transformed data as a list
  result = list(transformed_x = transformed_x, transformed_y = transformed_y)
  return(result)
}

#' Estimate the gamma rate that minimizes the AUC difference
#'
#' This function uses the bisection method to find the gamma distribution rate that minimizes the difference between
#' the computed AUC and the target AUC.
#'
#' @param y_population A numeric vector of observed values for the population Y.
#' @param rate_x The rate parameter of the gamma distribution for X.
#' @param auc_target The target AUC value to match.
#' @param tol The tolerance level for AUC difference (default is 0.01).
#' @param max_iter The maximum number of iterations for the bisection method (default is 100).
#' @param lower The lower bound used to start the algorithm.
#' @param upper The lower bound used to start the algorithm.
#' @return The estimated gamma rate or FALSE if the target AUC is not reached within the maximum iterations.
#' @export
get_gamma_rate = function(y_population, rate_x, auc_target, tol = 0.01, max_iter = 100, lower = 0, upper = 15) {

  iter = 0
  
  # Perform bisection method to find the gamma rate
  while (iter < max_iter) {
    bisectriz = (lower + upper) / 2
    iter = iter + 1
    
    # Sample X from the gamma distribution
    x_population = rgamma(1e5, shape = bisectriz, rate = rate_x)
    
    # Compute the AUC for the current iteration
    auc_iter = as.numeric(pROC::auc(response = c(rep(1, 1e5), rep(0, 1e5)), predictor = c(y_population, x_population)))
    cat('AUC iteration ', iter, ': ', auc_iter, '\n')
    
    # Check if the AUC difference is within the tolerance
    if (abs(auc_iter - auc_target) < tol) {
      return(bisectriz)
    }
    
    # Update bounds based on the AUC comparison
    if (auc_iter > auc_target) {
      upper = bisectriz
    }
    if (auc_iter < auc_target) {
      lower = bisectriz
    }
  }
  
  # Return FALSE if the target AUC is not reached
  return(FALSE)
}


#' Estimate the mean parameter from a gaussian that minimizes the AUC difference
#'
#' This function uses the bisection method to find the gaussian distribution mean that minimizes the difference between
#' the computed AUC and the target AUC.
#'
#' @param y_population A numeric vector of observed values for the population Y.
#' @param stdev_x The sigma parameter for the population X.
#' @param auc_target The target AUC value to match.
#' @param tol The tolerance level for AUC difference (default is 0.01).
#' @param max_iter The maximum number of iterations for the bisection method (default is 100).
#' @param lower The lower bound used to start the algorithm.
#' @param upper The upper bound used to start the algorithm.
#' @return The estimated mean or FALSE if the target AUC is not reached within the maximum iterations.
#' @export
get_mux_bisection = function(y_population, stdev_x, auc_target, tol = 0.01, max_iter = 100, lower = 0, upper = 10) {

  iter = 0
  
  # Perform bisection method to find the gamma rate
  while (iter < max_iter) {
    bisectriz = (lower + upper) / 2
    iter = iter + 1
    
    # Sample X from the gaussian distribution
    x_population = rnorm(1e5, bisectriz, stdev_x)

    
    # Compute the AUC for the current iteration
    auc_iter = as.numeric(pROC::auc(response = c(rep(1, 1e5), rep(0, 1e5)), predictor = c(y_population, x_population)))
    cat('AUC iteration ', iter, ': ', auc_iter, '\n')
    
    # Check if the AUC difference is within the tolerance
    if (abs(auc_iter - auc_target) < tol) {
      return(bisectriz)
    }
    
    # Update bounds based on the AUC comparison
    if (auc_iter > auc_target) {
      upper = bisectriz
    }
    if (auc_iter < auc_target) {
      lower = bisectriz
    }
  }
  
  # Return FALSE if the target AUC is not reached
  return(FALSE)
}

calculate_auc_normal = function(controls, cases){
  nAUC=pnorm((mean(cases)-mean(controls))/(sd(controls)^2+sd(cases)^2)^0.5)
  return(nAUC)
}

calculate_auc_kernel = function(controls, cases, bandwidth_method = "optimal", 
                                kernel_function = "gaussian", box_cox = FALSE){
  
  # Box-Cox transformation if requested
  if (box_cox){
    transformed_data = apply_box_cox(controls, cases)
    controls = transformed_data$transformed_x
    cases = transformed_data$transformed_y
  }
  
  n_controls = length(controls)
  n_cases = length(cases)
  
  # Bandwidth selection methods
  if (bandwidth_method == 'optimal') {
    bandwidth_controls = 1.06 * sd(controls) * n_controls^(-1/5)
    bandwidth_cases = 1.06 * sd(cases) * n_cases^(-1/5)
  }
  else if (bandwidth_method == 'hscv') {
    bandwidth_controls = ks::hscv(controls)
    bandwidth_cases = ks::hscv(cases)
  }
  else if (bandwidth_method == 'iqr') {
    bandwidth_controls = 0.9 * min(sd(controls), IQR(controls) / 1.34) * (n_controls^(-1/5)) 
    bandwidth_cases = 0.9 * min(sd(cases), IQR(cases) / 1.34) * (n_cases^(-1/5))
  }
  else {
    stop("Invalid bandwidth method. Choose 'optimal', 'hscv', or 'iqr'")
  }
  
  # Combined bandwidth for kernel convolution
  h_combined = sqrt(bandwidth_controls^2 + bandwidth_cases^2)
  
  # Define kernel CDF functions
  if (kernel_function == "gaussian") {
    kernel_cdf = function(u) { pnorm(u) }
  }
  else if (kernel_function == "epanechnikov") {
    kernel_cdf = function(u) {
      result = (-u^3 + 15*u + 2*5^(3/2)) / (4*5^(3/2))
      result = result * (abs(u) <= sqrt(5)) + 1 * (u > sqrt(5))
      return(result)
    }
  }
  else {
    stop("Invalid kernel function. Choose 'gaussian' or 'epanechnikov'")
  }
  
  # Calculate AUC using the kernel formula
  auc = 0
  for (i in 1:n_controls) {
    for (j in 1:n_cases) {
      auc = auc + kernel_cdf((cases[j] - controls[i]) / h_combined)
    }
  }
  
  auc = auc / (n_controls * n_cases)

  return(auc)
}


calculate_youden_normal = function(controls, cases){
  cstar=((mean(cases)*sd(controls)^2-mean(controls)*sd(cases)^2)-sd(controls)*sd(cases)*sqrt((mean(controls)-mean(cases))^2+(sd(controls)^2-sd(cases)^2)*log(sd(controls)^2/sd(cases)^2)))/(sd(controls)^2-sd(cases)^2)
  youden=pnorm((cstar-mean(controls))/sd(controls))-pnorm((cstar-mean(cases))/sd(cases))
  return(youden)
}


calculate_youden_kernel = function(controls, cases, bandwidth_method, mesh_size_kernel = 1000, box_cox = FALSE){
    # Apply Box-Cox transformation if specified
  if (box_cox) {
    transformed_data = apply_box_cox(controls, cases)
    controls = transformed_data$transformed_x
    cases = transformed_data$transformed_y
  }
  
  # Combine the samples
  combined_sample = c(controls, cases)
  
  # Bandwidth selection methods
  if (bandwidth_method == 'optimal') {
    bandwidth_controls = 1.06 * sd(controls) * length(controls)^(-1/5)
    bandwidth_cases = 1.06 * sd(cases) * length(cases)^(-1/5)
  }
  else if (bandwidth_method == 'hscv') {
    bandwidth_controls = ks::hscv(controls)
    bandwidth_cases = ks::hscv(cases)
  }
  else if (bandwidth_method == 'iqr') {
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


  return(max(estimated_dist_cases - estimated_dist_controls))
  
}


