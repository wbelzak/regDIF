#' Regularized Differential Item Functioning
#'
#' Performs regularization of DIF effects in item response theory and confirmatory factor analysis models via penalized expectation-maximization.
#'
#' @usage
#' regDIF(x,
#'        y,
#'        family = c("bernoulli","categorical","gaussian"),
#'        penalty = c("lasso","mcp"),
#'        nlambda = 100,
#'        lambda.max = 2,
#'        gamma = 3,
#'        lambda = NULL,
#'        anchor = NULL,
#'        rasch = FALSE,
#'        standardize = TRUE,
#'        quadpts = 81,
#'        control = list())
#'
#' @param x Matrix or dataframe of DIF predictors.
#' @param y Matrix or dataframe of item responses. See below for supported distributions.
#' @param family Character value or vector indicating the item response distributions via \code{y}. For scales where item responses are of one type only, the user may input one character value indicating the type (e.g., \code{"categorical"}). For mixed item types, the user must specify a vector of characters in the order that corresponds to the response matrix via \code{y}; e.g., \code{c(rep("categorical",2)}\code{, "bernoulli"}\code{, rep("gaussian",3))}. Supports:
#' \itemize{
#'    \item{\code{"bernoulli"} - Bernoulli item response via logistic link function (i.e., 1PL or 2PL model, see rasch option below for 1PL). Must be numeric/integer (2 unique values), factor (2 levels), or logical.}
#'    \item{\code{"categorical"} - Categorical item response via ordered logistic link function (i.e., Graded Response Model). Must be numeric/integer or factor.}
#'    \item{\code{"gaussian"} - Gaussian item response via identity link function (i.e., Confirmatory Factor Analysis). Must be numeric/integer.}}
#' @param penalty Character value indicating the penalty function to use. Supports the least absolute selection and shrinkage operator (LASSO) and the minimax concave penalty (MCP).
#' @param nlambda Numeric value indicating how many lambda values to fit. Default is 100.
#' @param lambda.max Numberic value indicating the maximum lambda parameter to use for internal construction of lambda vector. Default is 3. Must be large enough to shrink all DIF effects to zero to begin with.
#' @param gamma Numeric value indicating the gamma parameter in the MCP function. Gamma controls the degree of tapering of DIF effects as lambda decreases. Larger gamma leads to faster tapering (less bias but possibly more unstable optimization), whereas smaller gamma leads to slower tapering (more bias but more stable optimization). Default is 3. Must be greater than 1.
#' @param lambda Optional numeric vector of lambda values \eqn{\ge} 0. If lambda is supplied, this overrides the automatic construction of lambda values via \code{nlambda}. Must be non-negative and in descending order, from largest to smallest values (e.g., \code{seq(1,0,-.01)}.
#' @param anchor Optional numeric value or vector indicating which item response(s) are anchors (e.g., \code{anchor = 1}). Default is \code{NULL}, meaning at least one DIF effect per covariate will be fixed to zero as lambda approaches 0 (required to identify the model).
#' @param rasch Logical value indicating whether to constrain item slopes to 1 (i.e., equal slopes). If \code{TRUE}, no slope DIF will be evaluated. Default is \code{FALSE}.
#' @param standardize Logical value indicating whether to standardize DIF covariates for regularization. Default is \code{TRUE}, as it is recommended that all covariates be on the same scale.
#' @param quadpts Numeric value indicating the number of quadrature points to be used in approximating the latent variable distribution during estimation. Default is \code{81}. Not recommended to be below 10 for accurate estimation.
#' @param control Optional list of optimization parameters. May be:
#' \describe{
#'    \item{tol}{Convergence threshold of EM algorithm. Default is \code{10^-5}.}
#'    \item{maxiter}{Maximum number of EM iterations. Default is \code{10000}.}}
#'
#' @return Function returns an object of class \code{regDIF}
#'
#' @examples
#' \dontrun{
#'
#' library(regDIF)
#' head(ida)
#' x <- ida[,7:9]
#' y <- ida[,1:6]
#' fit <- regDIF(x, y, family = "bernoulli", penalty = "lasso")
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
                   family = c("bernoulli","categorical","gaussian"),
                   penalty = c("lasso","mcp"),
                   nlambda = 100,
                   lambda.max = 2,
                   gamma = 3,
                   lambda = NULL,
                   anchor = NULL,
                   rasch = FALSE,
                   standardize = TRUE,
                   quadpts= 81,
                   control = list()){

  #preprocess data
  call <- match.call()
  data_scrub <- preprocess(x,y,family,penalty,nlambda,lambda.max,lambda,anchor,rasch,standardize,quadpts,control,call)

  #Run Reg-DIF by looping through lambda
  for(pen in 1:length(data_scrub$lambda)){

    #obtain regDIF estimates
    estimates <- em_estimation(data_scrub$p,data_scrub$responses,data_scrub$predictors,data_scrub$theta,data_scrub$itemtypes,penalty,data_scrub$lambda,gamma,pen,anchor,rasch,data_scrub$final.control,data_scrub$samp_size,data_scrub$num_items,data_scrub$num_responses,data_scrub$num_predictors,data_scrub$num_quadpts)

    #postprocess data
    data_final <- postprocess(estimates,data_scrub$responses,data_scrub$predictors,y,x,data_scrub$theta,data_scrub$lambda,pen,anchor,data_scrub$final.control,data_scrub$final,data_scrub$samp_size,data_scrub$num_responses,data_scrub$num_predictors,data_scrub$num_items,data_scrub$num_quadpts)

    #update parameter estimates for next lambda value
    data_scrub$p <- estimates[[2]]
    data_scrub$final <- data_final
  }

  #Obtain final results
  class(data_final) <- "regDIF"
  return(data_final)

}

