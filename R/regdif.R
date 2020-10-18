#' Regularized Differential Item Functioning
#'
#' Performs regularization of DIF effects in item response theory and confirmatory factor analysis models via penalized expectation-maximization.
#'
#' @usage
#' regDIF(item.data,
#'        predictor.data,
#'        item.type = c("bernoulli","categorical","gaussian"),
#'        penalty = c("lasso","mcp"),
#'        ntau = 100,
#'        tau.max = 2,
#'        alpha = 1,
#'        gamma = 3,
#'        tau = NULL,
#'        anchor = NULL,
#'        rasch = FALSE,
#'        impact.data = list(mean = NULL, var = NULL),
#'        standardize = TRUE,
#'        quadpts = 15,
#'        control = list())
#'
#' @param item.data Matrix or dataframe of item responses. See below for supported distributions.
#' @param predictor.data Matrix or dataframe of DIF and/or impact predictors. See below for option to specify different predictors for impact model.
#' @param item.type Character value or vector indicating the item response distributions via \code{y}. For scales where item responses are of one type only, the user may input one character value indicating the type (e.g., \code{"categorical"}). For mixed item types, the user must specify a vector of characters in the order that corresponds to the response matrix via \code{y}; e.g., \code{c(rep("categorical",2)}\code{, "bernoulli"}\code{, rep("gaussian",3))}. Supports:
#' \itemize{
#'    \item{\code{"bernoulli"} - Bernoulli item response via logistic link function (i.e., 1PL or 2PL model, see rasch option below for 1PL). Must be numeric/integer (2 unique values), factor (2 levels), or logical.}
#'    \item{\code{"categorical"} - Categorical item response via ordered logistic link function (i.e., Graded Response Model). Must be numeric/integer or factor.}
#'    \item{\code{"gaussian"} - Gaussian item response via identity link function (i.e., Confirmatory Factor Analysis). Must be numeric/integer.}}
#' @param penalty Character value indicating the penalty function to use. Supports:
#' \itemize{
#'    \item{\code{"lasso"} - The least absolute selection and shrinkage operator (LASSO) which controls DIF selection through \eqn{\tau} (tau).}
#'    \item{\code{"mcp"} - The minimax concave penalty (MCP) which controls DIF selection through \eqn{\tau} (tau) and estimator bias through \eqn{\gamma} (gamma).}}
#' @param ntau Numeric value indicating how many tau values to fit. Default is 100.
#' @param tau.max Numberic value indicating the maximum tau parameter to use for internal construction of tau vector. Default is 3. Must be large enough to shrink all DIF effects to zero to begin with.
#' @param alpha Numeric value indicating the alpha parameter in the elastic net penalty function. Alpha controls the degree to which LASSO or ridge is used during regularization. Default is 1, which is equivalent to LASSO. For ridge, set alpha to 0. NOTE: If using MCP penalty, alpha may not be exactly 0.
#' @param gamma Numeric value indicating the gamma parameter in the MCP function. Gamma controls the degree of tapering of DIF effects as tau decreases. Larger gamma leads to faster tapering (less bias but possibly more unstable optimization), whereas smaller gamma leads to slower tapering (more bias but more stable optimization). Default is 3. Must be greater than 1.
#' @param tau Optional numeric vector of tau values \eqn{\ge} 0. If tau is supplied, this overrides the automatic construction of tau values via \code{ntau}. Must be non-negative and in descending order, from largest to smallest values (e.g., \code{seq(1,0,-.01)}.
#' @param anchor Optional numeric value or vector indicating which item response(s) are anchors (e.g., \code{anchor = 1}). Default is \code{NULL}, meaning at least one DIF effect per covariate will be fixed to zero as tau approaches 0 (required to identify the model).
#' @param rasch Logical value indicating whether to constrain item slopes to 1 (i.e., equal slopes). If \code{TRUE}, no slope DIF will be evaluated. Default is \code{FALSE}.
#' @param impact.data Optional list of matrices or data frames with predictors for mean and variance impact. Allows for different sets of predictors on the mean and variance impact equations compared to the item response DIF equations.
#' @param standardize Logical value indicating whether to standardize DIF and impact covariates for regularization. Default is \code{TRUE}, as it is recommended that all covariates be on the same scale.
#' @param quadpts Numeric value indicating the number of quadrature points to be used in approximating the latent variable distribution during estimation. Uses adaptive quadrature. Default is \code{15}.
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
#' item.data <- ida[,1:6]
#' predictor.data <- ida[,7:9]
#' fit <- regDIF(item.data, predictor.data, item.type = "bernoulli", penalty = "lasso")
#' fit
#'
#' }
#'
#' @import stats utils
#' @importFrom Rcpp sourceCpp
#' @useDynLib regDIF, .registration = TRUE
#'
#' @export

regDIF <- function(item.data,
                   predictor.data,
                   item.type = c("bernoulli","categorical","gaussian"),
                   penalty = c("mcp","lasso"),
                   ntau = 100,
                   tau.max = 2,
                   alpha = 1,
                   gamma = 3,
                   tau = NULL,
                   anchor = NULL,
                   rasch = FALSE,
                   impact.data = list(mean = NULL, var = NULL),
                   standardize = TRUE,
                   quadpts = 15,
                   control = list()){

  # item.type <- "bernoulli";penalty="lasso";ntau=100;tau.max=2;alpha=1;gamma=3;tau=NULL;anchor=1;rasch=F;impact.data=list(mean = NULL, var = NULL);standardize=F;quadpts=15;control = list()

  #obtain larger tau if necessary
  need_larger_tau <- TRUE #only true to start
  tau_times <- 0
  while(need_larger_tau | tau_times < 6){

    #increase tau.max if all penalized parameters have NOT been removed from model (unless specifying anchor item)
    if(need_larger_tau == FALSE) break
    if(tau_times > 0) {
      tau.max <- tau.max*1.5
      tau[1] <- tau[1]*1.5
    }
    #if too many tau.max values have been tried, stop.
    if(tau_times == 5){
      print(coef(data_scrub$final))
      stop("tau.max is too small.\n  Three possible solutions:\n  1. Increase tau.max large enough to ensure all DIF parameters are removed from the model.\n  2. Standardize predictors if not already standardized.\n  3. Provide anchor item(s).", call. = TRUE)
    }

    #preprocess data
    call <- match.call()
    data_scrub <- preprocess(item.data,predictor.data,item.type,penalty,ntau,tau.max,tau,anchor,rasch,impact.data,standardize,quadpts,control,call)

    #Run Reg-DIF by looping through tau
    for(pen in 1:length(data_scrub$tau)){

      #obtain regDIF estimates
      estimates <- em_estimation(data_scrub$p,data_scrub$responses,data_scrub$predictors,data_scrub$mean_predictors,data_scrub$var_predictors,data_scrub$itemtypes,penalty,data_scrub$tau,alpha,gamma,pen,anchor,rasch,data_scrub$final.control,data_scrub$samp_size,data_scrub$num_items,data_scrub$num_responses,data_scrub$num_predictors,quadpts)

      #stop if tau.max is too small on first run
      p2 <- unlist(estimates[[1]])
      dif_parms <- p2[grep(paste0("cov"),names(p2))]
      if(any(data_scrub$itemtypes == "gaussian")) dif_parms <- dif_parms[-grep("s1",names(dif_parms))]
      if(is.null(anchor) & #no anchor
         pen == 1 & #first penalty value
         sum(abs(dif_parms)) > 0 & #not all DIF effects are zero
         alpha == 1 #alpha is 1 for lasso
         ){
        message("\nWarning: tau.max or user-defined tau value is too small to penalize all parameters to zero without anchor item. Automatically trying larger tau.max or tau value.")
        tau_times <- tau_times + 1
        break
      } else{
        need_larger_tau <- FALSE
      }

      #postprocess data
      data_final <- postprocess(estimates,data_scrub$responses,data_scrub$predictors,data_scrub$mean_predictors,data_scrub$var_predictors,item.data,predictor.data,impact.data,data_scrub$tau,alpha,pen,anchor,data_scrub$final.control,data_scrub$final,data_scrub$samp_size,data_scrub$num_responses,data_scrub$num_predictors,data_scrub$num_items,data_scrub$num_quadpts)

      #update parameter estimates for next tau value
      data_scrub$p <- estimates[[1]]
      data_scrub$final <- data_final
    }

  }

  #Obtain final results
  class(data_final) <- "regDIF"
  cat("\n")
  return(data_final)

}

