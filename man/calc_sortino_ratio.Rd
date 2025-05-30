% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/quant.R
\name{calc_sortino_ratio}
\alias{calc_sortino_ratio}
\title{Calculate Sortino Ratio}
\usage{
calc_sortino_ratio(x, rf, freq, mar = 0)
}
\arguments{
\item{x}{xts object}

\item{rf}{xts with risk-free time-series}

\item{freq}{days, weeks, months, quarters to annaulize return and risk}

\item{mar}{default `0`, sets min acceptable return to subset downside
volatility}
}
\value{
vector of Sortino Ratios
}
\description{
Calculate Sortino Ratio
}
\note{
The function doesn't test to see if `x` and `rf` are for intersecting
time freqs.
}
