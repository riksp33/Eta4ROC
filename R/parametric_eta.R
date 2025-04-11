#' Parametric Eta calculation from ROC curves
#'
#' This function calculates the parametric Eta based on the ROC curves of the control and case datasets. Optionally, it can apply a Box-Cox transformation to the data before performing the calculations.
#'
#' @param controls A numeric vector of control group data (e.g., healthy individuals).
#' @param cases A numeric vector of case group data (e.g., individuals with a condition).
#' @param box_cox A logical value indicating whether to apply a Box-Cox transformation to the data before calculation. Defaults to `FALSE`.
#' @param mesh A numeric vector defining the mesh grid used for interpolation. Defaults to `seq(0.00001, 0.99999, length.out = 10000)`, which creates a fine grid.
#'
#' @details 
#' The function computes the parametric Eta by calculating the ROC and its derivative for both the control and case datasets. If the `box_cox` parameter is set to `TRUE`, a Box-Cox transformation is applied to the control and case data prior to the Eta calculation. 
#' The Eta is derived using the relationship between the ROC curve and its derivative, which is then passed to the `eta_from_roc_curves` function.
#'
#' @return A numeric value representing the parametric Eta based on the ROC curves.
#' 
#' @examples
#' # Example usage without Box-Cox transformation
#' controls = rnorm(1000, mean = 0, sd = 1)
#' cases = rnorm(1000, mean = 1, sd = 1)
#' eta = parametric_eta(controls, cases)
#'
#' # Example usage with Box-Cox transformation
#' eta_boxcox = parametric_eta(controls, cases, box_cox = TRUE)
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