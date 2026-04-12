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

#' Partial AUC under normality — exact closed form
#'
#' ROC curve: TPR(p) = Phi(a + b * Phi^{-1}(p)),  a = (mu_y-mu_x)/sigma_y,  b = sigma_x/sigma_y
#' Closed form after substitution u = Phi^{-1}(p):
#'   pAUC(t0) = Phi_2( Phi^{-1}(t0), a/sqrt(1+b^2); rho = -b/sqrt(1+b^2) )
#'
#' @param controls Numeric vector of controls (or pass moments directly via mu_x, sigma_x).
#' @param cases    Numeric vector of cases.
#' @param t0       Upper FPR limit in (0, 1] (default 1 -> full AUC).
#' @return Exact partial AUC.
#' @export
calculate_pauc_normal = function(controls, cases, t0 = 1) {
  mu_x    <- mean(controls);  sigma_x <- sd(controls)
  mu_y    <- mean(cases);     sigma_y <- sd(cases)

  a   <- (mu_y - mu_x) / sigma_y
  b   <- sigma_x / sigma_y
  rho <- -b / sqrt(1 + b^2)

  mvtnorm::pmvnorm(
    lower = c(-Inf, -Inf),
    upper = c(qnorm(t0), a / sqrt(1 + b^2)),
    mean  = c(0, 0),
    sigma = matrix(c(1, rho, rho, 1), nrow = 2)
  ) |> as.numeric()
}


#' Partial Youden index under normality — exact closed form
#'
#' The unrestricted optimum c* solves a quadratic (existing calculate_youden_normal).
#' The partial optimum is max(c*, c_boundary) where c_boundary = qnorm(1-t0, mu_x, sigma_x).
#' If c* >= c_boundary the restriction is not binding and both coincide.
#' If c* <  c_boundary the Youden function is evaluated at c_boundary directly.
#'
#' @param controls Numeric vector of controls.
#' @param cases    Numeric vector of cases.
#' @param t0       Upper FPR limit in (0, 1] (default 1 -> full Youden).
#' @return Partial Youden index value.
#' @export
calculate_pyouden_normal = function(controls, cases, t0 = 1) {
  mu_x    <- mean(controls);  sigma_x <- sd(controls)
  mu_y    <- mean(cases);     sigma_y <- sd(cases)

  # Unrestricted optimal threshold (closed form from existing function)
  c_star <- ((mu_y * sigma_x^2 - mu_x * sigma_y^2) -
               sigma_x * sigma_y * sqrt(
                 (mu_x - mu_y)^2 +
                   (sigma_x^2 - sigma_y^2) * log(sigma_x^2 / sigma_y^2)
               )) / (sigma_x^2 - sigma_y^2)

  # Boundary threshold corresponding to FPR = t0
  c_boundary <- qnorm(1 - t0, mean = mu_x, sd = sigma_x)

  # Effective threshold: push to boundary if unrestricted optimum violates constraint
  c_eff <- max(c_star, c_boundary)

  # Youden at effective threshold
  tpr <- 1 - pnorm(c_eff, mean = mu_y, sd = sigma_y)
  fpr <- 1 - pnorm(c_eff, mean = mu_x, sd = sigma_x)

  return(tpr - fpr)
}

#' Partial AUC using kernel density estimation
#'
#' Restricts the AUC integral to the region where FPR <= t0. The threshold
#' c*(t0) is found via uniroot on the kernel CDF of controls, then only
#' controls above that threshold contribute to the sum.
#'
#' With t0 = 1 the result equals calculate_auc_kernel().
#'
#' @param controls         Numeric vector of controls.
#' @param cases            Numeric vector of cases.
#' @param t0               Upper FPR limit in (0, 1] (default 1 -> full AUC).
#' @param bandwidth_method One of 'optimal', 'hscv', 'iqr'.
#' @param kernel_function  'gaussian' or 'epanechnikov'.
#' @param box_cox          Logical; apply Box-Cox transformation.
#' @return Partial AUC value.
#' @export
calculate_pauc_kernel = function(controls, cases, t0 = 1,
                                  bandwidth_method = "optimal",
                                  kernel_function  = "gaussian",
                                  box_cox          = FALSE) {
  if (box_cox) {
    transformed = apply_box_cox(controls, cases)
    controls    = transformed$transformed_x
    cases       = transformed$transformed_y
  }

  n_controls = length(controls)
  n_cases    = length(cases)

  if (bandwidth_method == "optimal") {
    bw_x = 1.06 * sd(controls) * n_controls^(-1/5)
    bw_y = 1.06 * sd(cases)    * n_cases^(-1/5)
  } else if (bandwidth_method == "hscv") {
    bw_x = ks::hscv(controls)
    bw_y = ks::hscv(cases)
  } else if (bandwidth_method == "iqr") {
    bw_x = 0.9 * min(sd(controls), IQR(controls) / 1.34) * n_controls^(-1/5)
    bw_y = 0.9 * min(sd(cases),    IQR(cases)    / 1.34) * n_cases^(-1/5)
  } else {
    stop("Invalid bandwidth_method. Choose 'optimal', 'hscv', or 'iqr'.")
  }

  h_combined = sqrt(bw_x^2 + bw_y^2)

  if (kernel_function == "gaussian") {
    kernel_cdf = function(u) pnorm(u)
  } else if (kernel_function == "epanechnikov") {
    kernel_cdf = function(u) {
      result = (-u^3 + 15 * u + 2 * 5^(3/2)) / (4 * 5^(3/2))
      result * (abs(u) <= sqrt(5)) + 1 * (u > sqrt(5))
    }
  } else {
    stop("Invalid kernel_function. Choose 'gaussian' or 'epanechnikov'.")
  }

  # Kernel CDF of controls at a given threshold
  F_controls_hat = function(c_val) {
    mean(kernel_cdf((c_val - controls) / bw_x))
  }

  # Find c*(t0) such that FPR_hat(c*) = t0, i.e. F_controls_hat(c*) = 1 - t0
  if (t0 >= 1) {
    c_upper = -Inf
  } else {
    lo      = min(controls) - 5 * bw_x
    hi      = max(controls) + 5 * bw_x
    c_upper = uniroot(function(c) F_controls_hat(c) - (1 - t0),
                      lower = lo, upper = hi)$root
  }

  # Sum only over controls in the FPR <= t0 region (x >= c_upper)
  pauc = 0
  for (i in seq_len(n_controls)) {
    if (is.finite(c_upper) && controls[i] < c_upper) next
    for (j in seq_len(n_cases)) {
      pauc = pauc + kernel_cdf((cases[j] - controls[i]) / h_combined)
    }
  }

  pauc = pauc / (n_controls * n_cases)
  return(pauc)
}

#' Partial Youden index using kernel density estimation
#'
#' Restricts the optimisation to mesh points where FPR <= t0.
#'
#' With t0 = 1 the result equals calculate_youden_kernel().
#'
#' @param controls         Numeric vector of controls.
#' @param cases            Numeric vector of cases.
#' @param t0               Upper FPR limit in (0, 1] (default 1 -> full Youden).
#' @param bandwidth_method One of 'optimal', 'hscv', 'iqr'.
#' @param mesh_size_kernel Number of mesh points (default 1000).
#' @param box_cox          Logical; apply Box-Cox transformation.
#' @return Partial Youden index value.
#' @export
calculate_pyouden_kernel = function(controls, cases, t0 = 1,
                                     bandwidth_method  = "optimal",
                                     mesh_size_kernel  = 1000,
                                     box_cox           = FALSE) {
  if (box_cox) {
    transformed = apply_box_cox(controls, cases)
    controls    = transformed$transformed_x
    cases       = transformed$transformed_y
  }

  if (bandwidth_method == "optimal") {
    bw_x = 1.06 * sd(controls) * length(controls)^(-1/5)
    bw_y = 1.06 * sd(cases)    * length(cases)^(-1/5)
  } else if (bandwidth_method == "hscv") {
    bw_x = ks::hscv(controls)
    bw_y = ks::hscv(cases)
  } else if (bandwidth_method == "iqr") {
    bw_x = 0.9 * min(sd(controls), IQR(controls) / 1.34) * length(controls)^(-1/5)
    bw_y = 0.9 * min(sd(cases),    IQR(cases)    / 1.34) * length(cases)^(-1/5)
  } else {
    stop("Invalid bandwidth_method. Choose 'optimal', 'hscv', or 'iqr'.")
  }

  mesh = seq(min(c(controls, cases)),
             max(c(controls, cases)),
             length.out = mesh_size_kernel)

  dist_controls = kernel_distribution_estimation(sort(controls), mesh, bw_x)
  dist_cases    = kernel_distribution_estimation(sort(cases),    mesh, bw_y)

  fpr_mesh  = 1 - dist_controls
  tpr_mesh  = 1 - dist_cases
  in_region = fpr_mesh <= t0

  if (!any(in_region)) {
    warning("No mesh points fall within FPR <= t0. Returning NA.")
    return(NA_real_)
  }

  return(max(tpr_mesh[in_region] - fpr_mesh[in_region]))
}