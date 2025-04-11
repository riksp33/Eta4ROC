#' Parametric Eta calculation from ROC curves
#'
#' This function calculates the parametric Eta based on the ROC curves of the control and case datasets. Optionally, it can apply a Box-Cox transformation to the data before performing the calculations.
#'
#' @param controls A numeric vector of control group data (e.g., healthy individuals).
#' @param cases A numeric vector of case group data (e.g., individuals with a condition).
#' @param t0 A numeric value indicating the cutoff point for integration. Only mesh points up to t0 are considered. Defaults to 1.
#' @param box_cox A logical value indicating whether to apply a Box-Cox transformation to the data before calculation. Defaults to `FALSE`.
#' @param mesh A numeric vector defining the mesh grid used for interpolation. Defaults to `seq(0.00001, 0.99999, length.out = 10000)`, which creates a fine grid.
#'
#' @details
#' The function computes the parametric Eta by:
#' 1. Optionally applying a Box-Cox transformation to normalize the data
#' 2. Calculating means and standard deviations of both groups
#' 3. Computing the binormal ROC curve and its derivative
#' 4. Passing these values to the `eta_from_roc_curves` function
#' 
#' The calculation assumes that the data (potentially after transformation) follows a normal distribution. 
#' The ratio of standard deviations (ro) and standardized difference (delta) are key parameters in the computation.
#'
#' @return A numeric value representing the parametric Eta based on the ROC curves.
#'
#' @examples
#' \dontrun{
#' # Generate sample data with moderate separation
#' controls <- rnorm(1000, mean = 0, sd = 1)
#' cases <- rnorm(1000, mean = 1, sd = 1)
#' 
#' # Calculate eta without transformation
#' eta <- parametric_eta(controls, cases)
#' print(eta)
#' 
#' # Calculate eta with Box-Cox transformation (useful for skewed data)
#' eta_boxcox <- parametric_eta(controls, cases, box_cox = TRUE)
#' print(eta_boxcox)
#' 
#' # Using a different cutoff point
#' eta_t0_0.8 <- parametric_eta(controls, cases, t0 = 0.8)
#' print(eta_t0_0.8)
#' 
#' # Using a coarser mesh for faster computation
#' coarse_mesh <- seq(0.001, 0.999, length.out = 1000)
#' eta_coarse <- parametric_eta(controls, cases, mesh = coarse_mesh)
#' print(eta_coarse)
#' }
#'
#' @seealso \code{\link{eta_from_roc_curves}} for the underlying eta calculation, 
#'   \code{\link{apply_box_cox}} for details on the Box-Cox transformation
#'
#' @export
parametric_eta = function(controls, cases, t0 = 1, box_cox = FALSE, mesh = seq(0.00001, 0.99999, length.out = 10000)){

    if(box_cox){
        transformed_data = apply_box_cox(controls, cases)
        controls = transformed_data$transformed_x
        cases = transformed_data$transformed_y
    }

    mux=mean(controls)
    muy=mean(cases)

    sigmax=sd(controls)
    sigmay=sd(cases)

    ro=sigmax/sigmay
    delta=(muy-mux)/sigmay

    ROC=1-pnorm(qnorm(1-mesh,mux,sigmax),mux+delta*sigmax/ro,sigmax/ro)
    ROCprima=(ro*exp(-0.5*(delta+ro*qnorm(mesh))^2))/exp(-0.5*(qnorm(mesh))^2)

    return(eta_from_roc_curves(ROC, ROCprima, t0, mesh))
}