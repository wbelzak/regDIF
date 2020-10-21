#' Regularized Differential Item Functioning
#'
#' Identify DIF in item response theory models using regularization.
#'
#' @usage
#' regDIF(item.data,
#'        predictor.data,
#'        item.type = NULL,
#'        penalty.type = NULL,
#'        tau = NULL,
#'        num.tau = 100,
#'        max.tau = 2,
#'        alpha = 1,
#'        gamma = 3,
#'        anchor = NULL,
#'        impact.data = list(mean = NULL, var = NULL),
#'        standardize = TRUE,
#'        quadpts = 15,
#'        control = list())
#'
#' @param item.data Matrix or dataframe of item responses. See below for
#' supported distributions.
#' @param predictor.data Matrix or dataframe of DIF and/or impact predictors.
#' See below for option to specify different predictors for impact model.
#' @param item.type Optional character value or vector indicating the type of
#' item to be modeled. The default is NULL, corresponding to a 2PL or graded
#' item type. Different item types may be specified for a single model
#' by providing a vector equal in length to the number of items in item.data.
#' The options include:
#' \itemize{
#'    \item{\code{"rasch"} - Slopes constrained to 1 and intercepts freely
#'    estimated.}
#'    \item{\code{"2pl"} - Slopes and intercepts freely estimated.}
#'    \item{\code{"graded"} - Slopes, intercepts, and thresholds freely
#'    estimated.}
#'    \item{\code{"cfa"} - Not currently supported. Slopes and intercepts
#'    freely estimated with continuous (Gaussian) item responses.}}
#' @param penalty.type Optional character value indicating the penalty
#' function to use. The default is NULL, corresponding to the LASSO function.
#' The options include:
#' \itemize{
#'    \item{\code{"mcp"} - The minimax concave penalty (MCP), which controls
#'    DIF selection through \eqn{\tau} (tau) and estimator bias through
#'    \eqn{\gamma} (gamma).}
#'    \item{\code{"lasso"} - The least absolute selection and shrinkage
#'    operator (LASSO), which controls DIF selection through \eqn{\tau} (tau).}
#'    \item{\code{"ridge"} - The ridge penalty, which shrinks DIF effects (but
#'    does not select DIF) through \eqn{\tau} (tau).}
#'    \item{\code{"elastic.net"} - The elastic net penalty, which selects DIF
#'    effects using a combination of lasso and ridge penalties through
#'    \eqn{\tau} (tau) and \eqn{\alpha} (alpha). Using this option }}
#' @param tau Optional numeric vector of tau values \eqn{\ge} 0. If tau is
#' supplied, this overrides the automatic construction of tau values via
#' \code{max.tau}. Must be non-negative and in descending order, from largest
#' to smallest values (e.g., \code{seq(1,0,-.01)}.
#' @param num.tau Numeric value indicating how many tau values to fit. The
#' default is 100.
#' @param max.tau Numberic value indicating the maximum tau parameter to use
#' for internal construction of tau vector. Default is 2. Must be large enough
#' to shrink all DIF effects to zero to begin with.
#' @param alpha Numeric value indicating the alpha parameter in the elastic net
#' penalty function. Alpha controls the degree to which LASSO or ridge is used
#' during regularization. The default is 1, which is equivalent to LASSO.
#' NOTE: If using MCP penalty, alpha may not be exactly 0.
#' @param gamma Numeric value indicating the gamma parameter in the MCP
#' function. Gamma controls the degree of tapering of DIF effects as tau
#' decreases. Larger gamma leads to faster tapering (less bias but possibly
#' more unstable optimization), whereas smaller gamma leads to slower tapering
#' (more bias but more stable optimization). Default is 3. Must be greater
#' than 1.
#' @param anchor Optional numeric value or vector indicating which item
#' response(s) are anchors (e.g., \code{anchor = 1}). Default is \code{NULL},
#' meaning at least one DIF effect per covariate will be fixed to zero as
#' tau approaches 0 (required to identify the model).
#' @param impact.data Optional list of matrices or data frames with predictors
#' for mean and variance impact. Allows for different sets of predictors on
#' the mean and variance impact equations compared to the item response DIF
#' equations.
#' @param standardize Logical value indicating whether to standardize DIF and
#' impact covariates for regularization. Default is \code{TRUE}, as it is
#' recommended that all covariates be on the same scale.
#' @param quadpts Numeric value indicating the number of quadrature points to
#' be used in approximating the latent variable distribution during estimation.
#' Uses adaptive quadrature. Default is \code{15}.
#' @param control Optional list of optimization parameters. May be:
#' \describe{
#'    \item{tol}{Convergence threshold of EM algorithm. Default is
#'    \code{10^-5}.}
#'    \item{maxiter}{Maximum number of EM iterations. Default is
#'    \code{10000}.}}
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
#' fit <- regDIF(item.data, predictor.data)
#' summary(fit)
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
                   item.type = NULL,
                   penalty.type = NULL,
                   tau = NULL,
                   num.tau = 100,
                   max.tau = 2,
                   alpha = 1,
                   gamma = 3,
                   anchor = NULL,
                   impact.data = list(mean = NULL, var = NULL),
                   standardize = TRUE,
                   quadpts = 15,
                   control = list()) {


  # Obtain larger tau if necessary.
  need_larger_tau <- TRUE
  tau_times <- 0
  while(need_larger_tau | tau_times < 6){

    # Increase tau.max if all penalized parameters have NOT been removed from
    # model (unless specifying anchor item).
    if(need_larger_tau == FALSE) break
    if(tau_times > 0) {
      max.tau <- max.tau*1.5
      tau[1] <- tau[1]*1.5
    }
    # If too many tau.max values have been tried, stop.
    if(tau_times == 5){
      print(coef(data_scrub$final))
      stop("max.tau is too small.\n  Three possible solutions:\n  1. Increase
           max.tau large enough to ensure all DIF parameters are removed from
           the model.\n  2. Standardize predictors if not already
           standardized.\n  3. Provide anchor item(s).", call. = TRUE)
    }

    # Pre-process data.
    call <- match.call()
    data_scrub <- preprocess(item.data,
                             predictor.data,
                             item.type,
                             num.tau,
                             max.tau,
                             tau,
                             anchor,
                             rasch,
                             impact.data,
                             standardize,
                             quadpts,
                             control,
                             call)
    if(is.null(penalty.type)) penalty.type <- "lasso"

    # Run Reg-DIF by looping through tau.
    for(pen in 1:length(data_scrub$tau)){

      # Obtain regDIF estimates.
      estimates <- em_estimation(data_scrub$p,
                                 data_scrub$item.data,
                                 data_scrub$predictor.data,
                                 data_scrub$mean_predictors,
                                 data_scrub$var_predictors,
                                 data_scrub$item.type,
                                 penalty.type,
                                 data_scrub$tau,
                                 alpha,
                                 gamma,
                                 pen,
                                 anchor,
                                 data_scrub$final.control,
                                 data_scrub$samp_size,
                                 data_scrub$num_items,
                                 data_scrub$num_responses,
                                 data_scrub$num_predictors,
                                 quadpts)

      # Stop if tau.max is too small on first run.
      p2 <- unlist(estimates[[1]])
      dif_parms <- p2[grep(paste0("cov"),names(p2))]

      if(any(data_scrub$itemtypes == "continuous")) {
        dif_parms <- dif_parms[-grep("s1",names(dif_parms))]
      }


      if(is.null(anchor) &
         pen == 1 &
         sum(abs(dif_parms)) > 0 &
         alpha == 1
         ) {
        message("\nWarning: tau.max or user-defined tau value is too small
                to penalize all parameters to zero without anchor item.
                Automatically trying larger tau.max or tau value.")
        tau_times <- tau_times + 1
        break
      } else{
        need_larger_tau <- FALSE
      }

      # Post-process data.
      data_final <- postprocess(estimates,
                                data_scrub$item.data,
                                data_scrub$predictor.data,
                                data_scrub$mean_predictors,
                                data_scrub$var_predictors,
                                impact.data,
                                data_scrub$tau,
                                alpha,
                                pen,
                                anchor,
                                data_scrub$final.control,
                                data_scrub$final,
                                data_scrub$samp_size,
                                data_scrub$num_responses,
                                data_scrub$num_predictors,
                                data_scrub$num_items,
                                data_scrub$num_quadpts)

      # Update parameter estimates for next tau value.
      data_scrub$p <- estimates[[1]]
      data_scrub$final <- data_final
    }

  }

  # Obtain final results.
  class(data_final) <- "regDIF"
  cat("\n")
  return(data_final)

}

