#' Regularized Differential Item Functioning
#'
#' Identify DIF in item response theory models using regularization.
#'
#' @usage
#' regDIF(item.data,
#'        pred.data,
#'        item.type = NULL,
#'        pen.type = NULL,
#'        tau = NULL,
#'        num.tau = 100,
#'        alpha = 1,
#'        gamma = 3,
#'        anchor = NULL,
#'        stdz = TRUE,
#'        adapt.quad = FALSE,
#'        num.quad = if(adapt.quad) 15,
#'        control = list())
#'
#' @param item.data Matrix or data frame of item responses. See below for
#' supported item types.
#' @param pred.data Matrix or data frame of DIF/impact predictors.
#' See \code{control} option below to specify different predictors for impact
#' model.
#' @param item.type Optional character value or vector indicating the type of
#' item to be modeled. The default is NULL, corresponding to a 2PL or graded
#' item type. Different item types may be specified for a single model
#' by providing a vector equal in length to the number of items in item.data.
#' The options include:
#' \itemize{
#'    \item{\code{"Rasch"} - Slopes constrained to 1 and intercepts freely
#'    estimated.}
#'    \item{\code{"2PL"} - Slopes and intercepts freely estimated.}
#'    \item{\code{"Graded"} - Slopes, intercepts, and thresholds freely
#'    estimated.}}
#' @param pen.type Optional character value indicating the penalty
#' function to use. The default is NULL, corresponding to the LASSO function.
#' The options include:
#' \itemize{
#'    \item{\code{"lasso"} - The least absolute selection and shrinkage
#'    operator (LASSO), which controls DIF selection through \eqn{\tau} (tau).}
#'    \item{\code{"mcp"} - The minimax concave penalty (MCP), which controls
#'    DIF selection through \eqn{\tau} (tau) and estimator bias through
#'    \eqn{\gamma} (gamma). Uses the firm-thresholding penalty function.}}
#' @param tau Optional numeric vector of tau values \eqn{\ge} 0. If tau is
#' supplied, this overrides the automatic construction of tau values via
#' \code{max.tau}. Must be non-negative and in descending order, from largest
#' to smallest values (e.g., \code{seq(1,0,-.01)}.
#' @param num.tau Numeric value indicating how many tau values to fit. The
#' default is 100.
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
#' @param stdz Logical value indicating whether to standardize DIF and
#' impact predictors for regularization. Default is \code{TRUE}, as it is
#' recommended that all predictors be on the same scale.
#' @param adapt.quad Adaptive quadrature is not currently supported.
#' Logical value indicating whether to use adaptive quadrature. The default is
#' \code{FALSE}.
#' @param num.quad Numeric value indicating the number of quadrature points to
#' be used in approximating the latent variable distribution during estimation.
#' Default is \code{15} when using adaptive quadrature (not currently
#' supported). When using fixed-point quadrature, the default is
#' \code{21} points if all items are binary or \code{51} points if at least one
#' item is ordered categorical.
#' @param control Optional list of different model specifications and
#' optimization parameters. May be:
#' \describe{
#'    \item{impact.data}{List of matrices or data frames with predictors for
#'    mean and variance impact. Allows for different sets of predictors on the
#'    mean and variance impact equations compared to the item response DIF
#'    equations. Default is for all predictors to affect the mean and variance
#'    equations.}
#'    \item{tol}{Convergence threshold of EM algorithm. Default is
#'    \code{10^-5}.}
#'    \item{maxiter}{Maximum number of EM iterations. Default is
#'    \code{5000}.}}
#'
#' @return Function returns an object of class \code{regDIF}
#'
#' @examples
#' \dontrun{
#'
#' library(regDIF)
#' head(ida)
#' item.data <- ida[,1:6]
#' pred.data <- ida[,7:9]
#' fit <- regDIF(item.data, pred.data, num.tau = 10)
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
                   pred.data,
                   item.type = NULL,
                   pen.type = NULL,
                   tau = NULL,
                   num.tau = 100,
                   alpha = 1,
                   gamma = 3,
                   anchor = NULL,
                   stdz = TRUE,
                   adapt.quad = FALSE,
                   num.quad = if(adapt.quad) 15,
                   control = list()) {


    # Pre-process data.
    call <- match.call()
    data_scrub <- preprocess(item.data,
                             pred.data,
                             item.type,
                             num.tau,
                             tau,
                             anchor,
                             stdz,
                             num.quad,
                             control,
                             pen.type,
                             adapt.quad,
                             call)

    # Identify minimum tau that removes all DIF.
    if(data_scrub$id_tau) {
      max.tau.fit <- identify_tau(data_scrub$p,
                                  data_scrub$item.data,
                                  data_scrub$pred.data,
                                  data_scrub$mean_predictors,
                                  data_scrub$var_predictors,
                                  data_scrub$item.type,
                                  data_scrub$theta,
                                  data_scrub$pen.type,
                                  data_scrub$tau_vec,
                                  num.tau,
                                  alpha,
                                  gamma,
                                  anchor,
                                  data_scrub$final.control,
                                  data_scrub$samp_size,
                                  data_scrub$num_items,
                                  data_scrub$num_responses,
                                  data_scrub$num_predictors,
                                  data_scrub$num.quad,
                                  adapt.quad)

      p2 <- unlist(max.tau.fit$p)
      dif_parms <- p2[grep(paste0("cov"),names(p2))]


      if(is.null(anchor) & sum(abs(dif_parms)) > 0 & alpha == 1) {
        stop(paste0("\nWarning: Automatically-generated or user-defined tau",
                    "value is too small to penalize all parameters to zero",
                    "without anchor item. Larger values of tau are needed."))
      }

      # Update vector of tau values based on identification of minimum tau value
      # which removes all DIF from the model.
      data_scrub$tau_vec <- seq((max.tau.fit$max_tau)**(1/3),0,
                                length.out = data_scrub$num.tau)**3

      # Post-process data.
      data_final <- postprocess(max.tau.fit,
                                data_scrub$item.data,
                                data_scrub$pred.data,
                                data_scrub$mean_predictors,
                                data_scrub$var_predictors,
                                data_scrub$tau_vec,
                                alpha,
                                pen=1,
                                anchor,
                                control,
                                data_scrub$final.control,
                                data_scrub$final,
                                data_scrub$samp_size,
                                data_scrub$num_responses,
                                data_scrub$num_predictors,
                                data_scrub$num_items,
                                data_scrub$num.quad)

      # Update parameter estimates.
      data_scrub$p <- max.tau.fit$p
      data_scrub$final <- data_final
    }

    # Run Reg-DIF by looping through tau.
    start_regDIF <- ifelse(data_scrub$id_tau,2,1)
    for(pen in start_regDIF:data_scrub$num.tau){

      # Obtain regDIF estimates.
      estimates <- em_estimation(data_scrub$p,
                                 data_scrub$item.data,
                                 data_scrub$pred.data,
                                 data_scrub$mean_predictors,
                                 data_scrub$var_predictors,
                                 data_scrub$item.type,
                                 data_scrub$theta,
                                 data_scrub$pen.type,
                                 data_scrub$tau_vec,
                                 alpha,
                                 gamma,
                                 pen,
                                 anchor,
                                 data_scrub$final.control,
                                 data_scrub$samp_size,
                                 data_scrub$num_items,
                                 data_scrub$num_responses,
                                 data_scrub$num_predictors,
                                 data_scrub$num.quad,
                                 adapt.quad)

      if(!data_scrub$id_tau) {
      p2 <- unlist(estimates[[1]])
      dif_parms <- p2[grep(paste0("cov"),names(p2))]

      if(is.null(anchor) & pen == 1 & sum(abs(dif_parms)) > 0 & alpha == 1) {
        stop(paste0("\nWarning: Automatically-generated or user-defined tau",
                    "value is too small to penalize all parameters to zero",
                    "without anchor item. Larger values of tau are needed."))
      }
      }

      # Post-process data.
      data_final <- postprocess(estimates,
                                data_scrub$item.data,
                                data_scrub$pred.data,
                                data_scrub$mean_predictors,
                                data_scrub$var_predictors,
                                data_scrub$tau_vec,
                                alpha,
                                pen,
                                anchor,
                                control,
                                data_scrub$final.control,
                                data_scrub$final,
                                data_scrub$samp_size,
                                data_scrub$num_responses,
                                data_scrub$num_predictors,
                                data_scrub$num_items,
                                data_scrub$num.quad)

      # Update parameter estimates for next tau value.
      data_scrub$p <- estimates[[1]]
      data_scrub$final <- data_final
    }


  # Obtain final results.
  class(data_final) <- "regDIF"
  cat("\n")
  return(data_final)

}

