% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/quant.R
\name{calc_sharpe_ratio}
\alias{calc_sharpe_ratio}
\title{Calculate Sharpe Ratio}
\usage{
calc_sharpe_ratio(x, rf, freq)
}
\arguments{
\item{x}{xts object}

\item{rf}{xts with risk-free time-series}

\item{freq}{days, weeks, months, quarters to annaulize return and risk}
}
\value{
vector of Sharpe Ratios
}
\description{
Calculate Sharpe Ratio
}
\note{
The function doesn't test to see if `x` and `rf` are for intersecting
time freqs.
}
