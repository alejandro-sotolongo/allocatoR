% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util_time_series.R
\name{eo_month}
\alias{eo_month}
\title{Convert dates to end of month}
\usage{
eo_month(x)
}
\arguments{
\item{x}{date or vector of dates}
}
\value{
dates shifted to end of month
}
\description{
Convert dates to end of month
}
\examples{
eo_month(as.Date("2024-11-29"))
}
