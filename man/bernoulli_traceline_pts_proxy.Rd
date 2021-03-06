% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/traceline.r
\name{bernoulli_traceline_pts_proxy}
\alias{bernoulli_traceline_pts_proxy}
\title{Binary item tracelines for proxy scores.}
\usage{
bernoulli_traceline_pts_proxy(p_item, prox_data, pred_data)
}
\arguments{
\item{p_item}{Vector of item parameters.}

\item{prox_data}{Vector of observed proxy scores.}

\item{pred_data}{Matrix or dataframe of DIF and/or impact predictors.}
}
\value{
a \code{"matrix"} of probability values for Bernoulli item likelihood using observed
proxy scores
}
\description{
Binary item tracelines for proxy scores.
}
\keyword{internal}
