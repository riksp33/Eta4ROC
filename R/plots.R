#' Plot ROC Curve
#'
#' This function creates and plots a ROC (Receiver Operating Characteristic) curve based on the given predicted values and true labels.
#' It combines the two input vectors (predicted and true labels) and computes the ROC curve using the `pROC` package.
#'
#' @param X vector of cases group.
#' @param Y vector of controls group.
#' @return A plot of the ROC curve.
#' @importFrom pROC roc
#' @export
roc_plot = function(Y, X){

  # Combinar los casos y controles en un solo vector
  pred = c(X, Y)
  
  # Etiquetas: 1 para casos, 0 para controles
  clases = c(rep(1, length(X)), rep(0, length(Y)))
  
  # Crear la curva ROC
  roc_obj = pROC::roc(clases, pred)
  
  # Graficar
  plot(roc_obj, col = "blue", main = "Curva ROC")
}
