
#' Calculate eta value from ROC curves
#'
#' This function calculates the eta value using the given ROC curves and mesh. The eta is computed based on the formula
#' that involves the differences between the ROC curve values and mesh intervals.
#'
#' @param roc A numeric vector of ROC curve values.
#' @param roc_dx A numeric vector of ROC curve derivatives (roc_prima).
#' @param mesh A numeric vector representing the mesh/grid used for integration.
#' @return The standardized eta value.
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


#' General function to calculate eta analitically for gaussian and gamma distributions
#'
#' This function allows calculating eta using custom distributions and their respective transformations.
#'
#' @param mux A numeric value for the mean of the distribution (e.g., mean for Normal, Log-Normal).
#' @param mesh_size The number of intervals to use for the mesh/grid (default is m).
#' @param param_1_y Either mu or shape (gaussian / gamma) of the non-diseased population.
#' @param param_2_y Either stdev or rate (gaussian / gamma) of the non-diseased population.
#' @param param_1_x Either mu or shape (gaussian / gamma) of the diseased population.
#' @param param_2_x Either stdev or rate (gaussian / gamma) of the diseased population.

#' @return The calculated eta value.
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