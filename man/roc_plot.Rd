% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plots.R
\name{roc_plot}
\alias{roc_plot}
\title{Plot ROC Curve}
\usage{
roc_plot(Y, X)
}
\arguments{
\item{Y}{vector of controls group.}

\item{X}{vector of cases group.}
}
\value{
A plot of the ROC curve.
}
\description{
This function creates and plots a ROC (Receiver Operating Characteristic) curve based on the given predicted values and true labels.
It combines the two input vectors (predicted and true labels) and computes the ROC curve using the \code{pROC} package.
}
