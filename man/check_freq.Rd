% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util_time_series.R
\name{check_freq}
\alias{check_freq}
\title{Check Frequency Input}
\usage{
check_freq(freq)
}
\arguments{
\item{freq}{frequency string: days, weeks, months, quarters, or years}
}
\value{
if valid frequency, else throws error
}
\description{
Check Frequency Input
}
\examples{
check_freq("Months")
check_freq("M")
}
