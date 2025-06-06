% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/quant.R
\name{risk_wgt}
\alias{risk_wgt}
\title{Calculate Volatility or TE Weights}
\usage{
risk_wgt(x, cov_mat)
}
\arguments{
\item{x}{vector of weights corresponding to `cov_mat`}

\item{cov_mat}{covariance matrix}
}
\value{
a vector of risk weights
}
\description{
Calculate Volatility or TE Weights
}
\details{
For TE `x` is a vector of active weights (weight - benchmark weight) and
the `cov_mat` is a union of the portfolio and benchmark positions that
corresponds to the active weight vector `x`
}
