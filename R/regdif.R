#' Regularized Differential Item Functioning
#'
#' Regularization of DIF parameters in item response theory (IRT) and moderated nonlinear factor analysis (MNLFA) models.
#'
#' @param x Matrix or dataframe of predictors (i.e., DIF covariates). Supports categorical and continuous predictors.
#' @param y Matrix or dataframe of item responses. Supports Bernoulli (e.g., 0,1), categorical (e.g., 0,1,2,...,k), and Gaussian item responses.
#' @param itemtypes Character value or vector indicating the item response distributions. For scales where item responses are of one type only, the user may input one character value indicating the type (e.g., \code{"categorical"}). For mixed item types, the user must specify a vector of characters in the order that corresponds to the response matrix (e.g. \code{c(rep("categorical",3),rep("gaussian",3))}). Supports categorical and continuous responses, including Bernoulli with logistic link (i.e., binary outcomes), categorical with ordered logistic link (i.e., graded response model), and Gaussian with identity link (i.e., factor analysis).
#' @param penalty Character value indicating the penalty function to use. Supports the least absolute selection and shrinkage operator (lasso) and the minimax concave penalty (mcp).
#' @param nlambda Numeric value indicating the number of lambda values to fit.
#' @param lambda.max Numberic value indicating the maximum lambda (shrinkage tuning parameter) to fit. Default is 3.
#' @param gamma Numeric value indicating the gamma (tapering tuning parameter) to fit. Default is 3.
#' @param lambda Optional numeric vector of tuning parameter values \eqn{\ge} 0. If supplied, must be in descending order, from largest to smallest values (e.g., \code{seq(1,0,-.01)}.
#' @param anchor Optional numeric vector indicating which items are anchors (e.g., \code{anchor = 1}). Default is \code{NULL}, meaning at least one DIF parameter (per covariate) will be fixed to zero as lambda approaches 0. This is required to identify the model.
#' @param rasch Logical value indicating whether to constrain item slopes to 1 (i.e., equal slopes). If \code{TRUE}, no slope DIF will be evaluated. Default is \code{FALSE}.
#' @param control Optional list of control parameters. See documentation.
#'
#' @return Function returns an object of class \code{regDIF}
#'
#' @examples
#' \dontrun{
#'
#' library(regDIF)
#' head(ida)
#' y <- ida[,1:6]
#' x <- ida[,7:9]
#' fit <- regDIF(x, y, itemtypes = "bernoulli", penalty = "lasso")
#' fit
#'
#' }
#'
#' @import stats utils
#' @importFrom Rcpp sourceCpp
#' @useDynLib regDIF, .registration = TRUE
#'
#' @export

regDIF <- function(x,
                   y,
                   itemtypes = c("bernoulli","categorical","gaussian"),
                   penalty = c("lasso","mcp"),
                   nlambda = 100,
                   lambda.max = 3,
                   gamma = 3,
                   lambda = NULL,
                   anchor = NULL,
                   rasch = FALSE,
                   control = list()){


  #preprocess data
  call <- match.call()
  data_scrub <- preprocess(x,y,itemtypes,penalty,nlambda,lambda.max,lambda,anchor,rasch,control,call)

  #Run Reg-DIF by looping through lambda
  for(pen in 1:length(lambda)){

    #obtain regDIF estimates
    estimates <- em_estimation(data_scrub$p,data_scrub$responses,data_scrub$predictors,data_scrub$theta,data_scrub$itemtypes,penalty,data_scrub$lambda,gamma,pen,anchor,rasch,data_scrub$final.control,data_scrub$samp_size,data_scrub$num_items,data_scrub$num_responses,data_scrub$num_predictors)

    data_final <- postprocess(estimates,data_scrub$responses,data_scrub$predictors,y,x,data_scrub$theta,lambda,pen,anchor,data_scrub$final.control,data_scrub$final,data_scrub$samp_size,data_scrub$num_responses,data_scrub$num_predictors,data_scrub$num_items)

  }

  #Obtain final results
  class(data_final) <- "regDIF"
  return(data_final)

}

