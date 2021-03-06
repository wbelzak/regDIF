% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/derivatives.r
\name{d_bernoulli_itemblock}
\alias{d_bernoulli_itemblock}
\title{Partial derivatives for binary items by item-blocks.}
\usage{
d_bernoulli_itemblock(
  p_item,
  etable,
  theta,
  pred_data,
  item_data_current,
  samp_size,
  num_items,
  num_predictors,
  num_quad
)
}
\arguments{
\item{p_item}{Vector of item parameters.}

\item{etable}{E-table for item.}

\item{theta}{Matrix of adaptive theta values.}

\item{pred_data}{Matrix or dataframe of DIF and/or impact predictors.}

\item{item_data_current}{Vector of current item responses.}

\item{samp_size}{Sample size in dataset.}

\item{num_items}{Number of items in dataset.}

\item{num_predictors}{Number of predictors in dataset.}

\item{num_quad}{Number of quadrature points used for approximating the
latent variable.}
}
\value{
a \code{"list"} of first and second partial derivatives for Bernoulli item likelihood (to
use with multivariate Newton-Raphson)
}
\description{
Partial derivatives for binary items by item-blocks.
}
\keyword{internal}
