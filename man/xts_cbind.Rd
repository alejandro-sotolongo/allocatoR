% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util_time_series.R
\name{xts_cbind}
\alias{xts_cbind}
\title{Column bind xts objects while preserving columns names}
\usage{
xts_cbind(x, y)
}
\arguments{
\item{x}{xts object}

\item{y}{xts object}
}
\value{
xts with `cbind(x, y)` with original column names of `x` and `y`
}
\description{
Column bind xts objects while preserving columns names
}
\details{
Column names will get converted to `data.frame` formats when `cbind` is
called on the xts object. E.g., Small Cap will be Small.Cap. This method
preserves the spaces or special characters in the original column names.
}
