% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util_time_series.R
\name{xts_to_tidy}
\alias{xts_to_tidy}
\title{Convert xts to tidy (long) data.frame}
\usage{
xts_to_tidy(x)
}
\arguments{
\item{x}{xts}
}
\description{
Convert xts to tidy (long) data.frame
}
\examples{
data(assets)
tidy_df <- xts_to_tidy(x)
head(tidy_df)
}
