#' Penalized expectation-maximization algorithm.
#'
#' @param p List of parameters with starting values obtained from preprocess.
#' @param item.data Matrix or dataframe of item responses.
#' @param predictor.data Matrix or dataframe of DIF and/or impact predictors.
#' @param mean_predictors Possibly different matrix of predictors for the mean
#' impact equation.
#' @param var_predictors Possibly different matrix of predictors for the
#' variance impact equation.
#' @param item.type Optional character value or vector indicating the type of
#' item to be modeled.
#' @param penalty.type Character value indicating the penalty function to use.
#' @param ntau Numeric value of how many to tau values to fit.
#' @param tau.max Numberic value indicating the maximum tau parameter.
#' @param alpha Numeric value indicating the alpha parameter in the elastic net
#' penalty function.
#' @param gamma Numeric value indicating the gamma parameter in the MCP
#' function.
#' @param tau Optional numeric vector of tau values.
#' @param anchor Optional numeric value or vector indicating which item
#' response(s) are anchors (e.g., \code{anchor = 1}).
#' @param rasch Logical value indicating whether to constrain item slopes
#' to 1 (i.e., equal slopes).
#' @param impact.data Optional list of matrices or data frames with predictors
#' for mean and variance impact.
#' @param standardize Logical value indicating whether to standardize DIF and
#' impact covariates for regularization.
#' @param quadpts Numeric value indicating the number of quadrature points.
#' @param control Optional list of optimization parameters.
#' @param call Defined from regDIF.
#'
#' @keywords internal
#'
em_estimation <- function(p,
                          item.data,
                          predictor.data,
                          mean_predictors,
                          var_predictors,
                          item.type,
                          penalty.type,
                          tau,
                          alpha,
                          gamma,
                          pen,
                          anchor,
                          final.control,
                          samp_size,
                          num_items,
                          num_responses,
                          num_predictors,
                          quadpts) {

  # Maximization settings.
  lastp <- p
  eps <- Inf
  iter = 0
  num_quadpts <- quadpts

  # Loop until convergence or maximum number of iterations.
  while(eps > final.control$tol & iter < final.control$maxit){

    # E-step: Evaluate Q function with current parameter estimates p.
    elist <- Estep(p,
                   item.data,
                   predictor.data,
                   mean_predictors,
                   var_predictors,
                   samp_size,
                   num_items,
                   num_responses,
                   num_quadpts)

    # M-step: Optimize parameters.
    p <- Mstep(p,
               item.data,
               predictor.data,
               mean_predictors,
               var_predictors,
               elist,
               item.type,
               penalty.type,
               tau[pen],
               alpha,
               gamma,
               anchor,
               samp_size,
               num_responses,
               num_items,
               num_quadpts,
               num_predictors)

    # Update and check for convergence: Calculate the difference
    # in parameter estimates from current to previous.
    eps = sqrt(sum((unlist(p)-unlist(lastp))^2))

    # Update parameter list.
    lastp <- p

    # Update the iteration number.
    iter = iter + 1
    if(iter == final.control$maxit) warning("EM iteration limit reached without convergence")

    cat('\r',sprintf("Models Completed: %d of %d  Iteration: %d  Change: %f", pen-1, length(tau), iter, round(eps,6))) #print information about optimization

    utils::flush.console()

  }

  # Get information criteria.
  infocrit <- information_criteria(elist,
                                   p,
                                   item.data,
                                   predictor.data,
                                   mean_predictors,
                                   var_predictors,
                                   elist$theta,
                                   tau[pen],
                                   gamma,
                                   samp_size,
                                   num_responses,
                                   num_items,
                                   num_quadpts)

  return(list(p,infocrit))

}
