#' Penalized expectation-maximization algorithm.
#'
#' @param p List of parameters with starting values obtained from preprocess.
#' @param item.data Matrix or data frame of item responses.
#' @param pred.data Matrix or data frame of DIF and/or impact predictors.
#' @param mean_predictors Possibly different matrix of predictors for the mean
#' impact equation.
#' @param var_predictors Possibly different matrix of predictors for the
#' variance impact equation.
#' @param item.type Optional character value or vector indicating the type of
#' item to be modeled.
#' @param theta Vector of fixed quadrature points.
#' @param pen.type Character value indicating the penalty function to use.
#' @param tau_vec Vector of tau values that either are automatically generated
#' or provided by the user. The first \code{tau_vec} will be equal to \code{Inf}
#' to identify a minimal value of tau in which all DIF is removed from the
#' model.
#' @param alpha Numeric value indicating the alpha parameter in the elastic net
#' penalty function.
#' @param gamma Numeric value indicating the gamma parameter in the MCP
#' function.
#' @param pen Index for the tau vector.
#' @param anchor Optional numeric value or vector indicating which item
#' response(s) are anchors (e.g., \code{anchor = 1}).
#' @param final.control Control parameters.
#' @param samp_size Numeric value indicating the sample size.
#' @param num_items Numeric value indicating the number of items.
#' @param num_responses Vector with number of responses for each item.
#' @param num_predictors Numeric value indicating the number of predictors.
#' @param num.quad Numeric value indicating the number of quadrature points.
#' @param adapt.quad Logical value indicating whether to use adaptive quad.
#' needs to be identified.
#'
#' @keywords internal
#'
em_estimation <- function(p,
                          item.data,
                          pred.data,
                          mean_predictors,
                          var_predictors,
                          item.type,
                          theta,
                          pen.type,
                          tau_vec,
                          alpha,
                          gamma,
                          pen,
                          anchor,
                          final.control,
                          samp_size,
                          num_items,
                          num_responses,
                          num_predictors,
                          num.quad,
                          adapt.quad) {

  # Maximization settings.
  lastp <- p
  eps <- Inf
  iter <- 0

  # Loop until convergence or maximum number of iterations.
  while(eps > final.control$tol && iter < final.control$maxit){


    # E-step: Evaluate Q function with current parameter estimates p.
    etable <- Estep(p,
                   item.data,
                   pred.data,
                   mean_predictors,
                   var_predictors,
                   theta,
                   samp_size,
                   num_items,
                   num_responses,
                   adapt.quad,
                   num.quad)

    # M-step: Optimize parameters.
    p <- Mstep(p,
               item.data,
               pred.data,
               mean_predictors,
               var_predictors,
               etable,
               item.type,
               pen.type,
               tau_vec[pen],
               alpha,
               gamma,
               anchor,
               samp_size,
               num_responses,
               num_items,
               num.quad,
               num_predictors)

    # Update and check for convergence: Calculate the difference
    # in parameter estimates from current to previous.
    eps = sqrt(sum((unlist(p)-unlist(lastp))^2))

    # Update parameter list.
    lastp <- p

    # Update the iteration number.
    iter = iter + 1
    if(iter == final.control$maxit) {
      warning("EM iteration limit reached without convergence")
    }

    cat('\r', sprintf("Models Completed: %d of %d  Iteration: %d  Change: %f",
                     pen-1,
                     length(tau_vec),
                     iter,
                     round(eps, nchar(final.control$tol))))

    utils::flush.console()


  }



  # Get information criteria.
  infocrit <- information_criteria(etable,
                                   p,
                                   item.data,
                                   pred.data,
                                   mean_predictors,
                                   var_predictors,
                                   gamma,
                                   samp_size,
                                   num_responses,
                                   num_items,
                                   num.quad)

  return(list(p,infocrit))

}
