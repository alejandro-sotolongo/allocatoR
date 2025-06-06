% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util_time_series.R
\name{colProd}
\alias{colProd}
\title{Calculate product of each column}
\usage{
colProd(x, add_1 = TRUE, na_rm = TRUE)
}
\arguments{
\item{x}{matrix, dataframe, or xts with numeric columns}

\item{add_1}{boolean to add 1 for geometric product of returns}

\item{na_rm}{boolean to remove NA values, passed through to na.rm of prod
function}
}
\value{
vector of column products
}
\description{
Calculate product of each column
}
\examples{
data(assets)
colProd(x)
}
