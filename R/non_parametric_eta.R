#' Non-parametric Eta estimation
#'
#' Estimates Eta directly from the control and case samples using the empirical
#' distribution of cases between consecutive control observations.
#'
#' @param controls A numeric vector of control group observations.
#' @param cases A numeric vector of case group observations.
#' @param t0 A numeric cutoff between 0 and 1. The default is 1.
#' @return The standardized non-parametric Eta estimate.
#' @details The empirical Eta is calculated from the counts of case observations
#'   between consecutive ordered controls, using the first
#'   \code{floor(t0 * (length(controls) + 1))} counts from the upper tail.
#'   The result is standardized with the package's Eta standardization.
#' @examples
#' controls <- c(0, 0.5, 2.5)
#' cases <- c(1, 2, 3, 4)
#' non_parametric_eta(controls, cases)
#' non_parametric_eta(controls, cases, t0 = 0.5)
#' @export
non_parametric_eta = function(controls, cases, t0 = 1) {
nd=length(cases)
ndbar=length(controls)
controls=sort(controls)
N=numeric()
for (k in 2:ndbar){
N[k]=sum(cases<=controls[k])-sum(cases<=controls[k-1])
}
N[1]=sum(cases<=controls[1])
N[ndbar+1]=sum(cases>controls[ndbar])
lim=floor(t0*(ndbar+1))
if (lim==0){
eta=0
} else {
N=N[(ndbar+2-lim):(ndbar+1)]
eta=sum((N/nd-1/(ndbar+1))^2)
}
return(standarization(eta))
}
