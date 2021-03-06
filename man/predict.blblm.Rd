% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/blblm.R
\name{predict.blblm}
\alias{predict.blblm}
\title{Model predictions on our blblm}
\usage{
\method{predict}{blblm}(object, new_data, confidence = FALSE, level = 0.95, ...)
}
\arguments{
\item{object}{object / model we want to test}

\item{new_data}{New data that we want to test our regression on}

\item{confidence}{Set to \code{TRUE} if you want a confidence
interval prediction.}

\item{level}{Confidence level we want to investigate}

\item{...}{Accepts any other argument you might want to pass.}
}
\description{
Model predictions on our blblm
}
\references{
http://www.sthda.com/english/articles/40-regression-analysis/166-predict-in-r-model-predictions-and-confidence-intervals/
}
