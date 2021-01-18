#' Penalized expectation-maximization algorithm.
#'
#' @param p List of parameters with starting values obtained from preprocess.
#' @param item_data Matrix or data frame of item responses.
#' @param pred_data Matrix or data frame of DIF and/or impact predictors.
#' @param mean_predictors Possibly different matrix of predictors for the mean
#' impact equation.
#' @param var_predictors Possibly different matrix of predictors for the
#' variance impact equation.
#' @param item_type Optional character value or vector indicating the type of
#' item to be modeled.
#' @param theta Vector of fixed quadrature points.
#' @param pen_type Character value indicating the penalty function to use.
#' @param tau_vec Vector of tau values that either are automatically generated
#' or provided by the user. The first \code{tau_vec} will be equal to \code{Inf}
#' to identify a minimal value of tau in which all DIF is removed from the
#' model.
#' @param id_tau Logical indicating whether to identify the minimum value of tau
#' in which all DIF parameters are removed from the model.
#' @param alpha Numeric value indicating the alpha parameter in the elastic net
#' penalty function.
#' @param gamma Numeric value indicating the gamma parameter in the MCP
#' function.
#' @param pen Index for the tau vector.
#' @param anchor Optional numeric value or vector indicating which item
#' response(s) are anchors (e.g., \code{anchor = 1}).
#' @param final_control Control parameters.
#' @param samp_size Numeric value indicating the sample size.
#' @param num_items Numeric value indicating the number of items.
#' @param num_responses Vector with number of responses for each item.
#' @param num_predictors Numeric value indicating the number of predictors.
#' @param num_quad Numeric value indicating the number of quadrature points.
#' @param adapt_quad Logical value indicating whether to use adaptive quad.
#' needs to be identified.
#' @param optim_method Character value indicating the type of optimization
#' method to use.
#' @param em_history List to save EM iterations for supplemental EM algorithm.
#'
#' @keywords internal
#'
em_estimation <- function(p,
                          item_data,
                          pred_data,
                          mean_predictors,
                          var_predictors,
                          item_type,
                          theta,
                          pen_type,
                          tau_vec,
                          id_tau,
                          alpha,
                          gamma,
                          pen,
                          anchor,
                          final_control,
                          samp_size,
                          num_items,
                          num_responses,
                          num_predictors,
                          num_quad,
                          adapt_quad,
                          optim_method,
                          em_history) {

  # Maximization settings.
  lastp <- p
  eps <- Inf
  iter <- 1

  # Loop until convergence or maximum number of iterations.
  while(eps > final_control$tol && iter < final_control$maxit){


    # E-step: Evaluate Q function with current parameter estimates p.
    etable <- Estep(p,
                    item_data,
                    pred_data,
                    mean_predictors,
                    var_predictors,
                    theta,
                    samp_size,
                    num_items,
                    num_responses,
                    adapt_quad,
                    num_quad)

    if(optim_method == "multi") {
      # M-step: Optimize parameters using multivariate NR.
      p <- Mstep_block(p,
                       item_data,
                       pred_data,
                       mean_predictors,
                       var_predictors,
                       etable,
                       item_type,
                       pen_type,
                       tau_vec[pen],
                       alpha,
                       gamma,
                       anchor,
                       samp_size,
                       num_responses,
                       num_items,
                       num_quad,
                       num_predictors)
    } else if(optim_method == "uni") {
      # M-step: Optimize parameters using one round of coordinate descent.
      p <- Mstep_cd(p,
                    item_data,
                    pred_data,
                    mean_predictors,
                    var_predictors,
                    etable,
                    item_type,
                    pen_type,
                    tau_vec[pen],
                    alpha,
                    gamma,
                    anchor,
                    samp_size,
                    num_responses,
                    num_items,
                    num_quad,
                    num_predictors)
    }

    # Update and check for convergence: Calculate the difference
    # in parameter estimates from current to previous.
    eps = sqrt(sum((unlist(p)-unlist(lastp))^2))

    # Save parameter estimates for supplemental em.
    em_history[[pen]][,iter] <- unlist(p)
    if(eps > final_control$tol) {
      em_history[[pen]] <- cbind(em_history[[pen]],
                                 matrix(0,ncol=1,nrow=length(unlist(p))))
    }


    # Update parameter list.
    lastp <- p

    # Update the iteration number.
    iter = iter + 1
    if(iter == final_control$maxit) {
      warning("EM iteration limit reached without convergence")
    }

    cat('\r', sprintf("Models Completed: %d of %d  Iteration: %d  Change: %f",
                     pen-1,
                     length(tau_vec),
                     iter,
                     round(eps, nchar(final_control$tol))))

    utils::flush.console()


  }

  # Get information criteria.
  infocrit <- information_criteria(etable,
                                   p,
                                   item_data,
                                   pred_data,
                                   mean_predictors,
                                   var_predictors,
                                   gamma,
                                   samp_size,
                                   num_responses,
                                   num_items,
                                   num_quad)

  # Option to identify maximum value of tau which removes all DIF from model.
  if(id_tau) {

    # Final M-step.
    if(optim_method == "multi") {
      max_tau <- Mstep_block_idtau(p,
                                   item_data,
                                   pred_data,
                                   mean_predictors,
                                   var_predictors,
                                   etable,
                                   item_type,
                                   pen_type,
                                   tau_vec[1],
                                   alpha,
                                   gamma,
                                   anchor,
                                   samp_size,
                                   num_responses,
                                   num_items,
                                   num_quad,
                                   num_predictors)
    } else if(optim_method == "uni") {
      max_tau <- Mstep_cd_idtau(p,
                                item_data,
                                pred_data,
                                mean_predictors,
                                var_predictors,
                                etable,
                                item_type,
                                pen_type,
                                tau_vec[1],
                                alpha,
                                gamma,
                                anchor,
                                samp_size,
                                num_responses,
                                num_items,
                                num_quad,
                                num_predictors)
    }

    # Return model results for maximum tau value.
    return(list(p=p,infocrit=infocrit,max_tau=max_tau,em_history=em_history))

  } else {

    # Return model results for all other tau values.
    return(list(p=p,infocrit=infocrit,em_history=em_history))
  }

}
