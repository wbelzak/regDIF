#' Penalized expectation-maximization algorithm.
#'
#' @param p List of parameters with starting values obtained from preprocess.
#' @param item_data Matrix or data frame of item responses.
#' @param pred_data Matrix or data frame of DIF and/or impact predictors.
#' @param prox_data Vector of observed proxy scores.
#' @param mean_predictors Possibly different matrix of predictors for the mean
#' impact equation.
#' @param var_predictors Possibly different matrix of predictors for the
#' variance impact equation.
#' @param item_type Character value or vector indicating the type of
#' item to be modeled.
#' @param theta Vector of fixed quadrature points.
#' @param pen_type Character value indicating the penalty function to use.
#' @param tau_vec Vector of tau values that either are automatically generated
#' or provided by the user. The first \code{tau_vec} will be equal to \code{Inf}
#' to identify a minimal value of tau in which all DIF is removed from the
#' model.
#' @param id_tau Logical indicating whether to identify the minimum value of tau
#' in which all DIF parameters are removed from the model.
#' @param num_tau Numeric value indicating the number of tau values to run
#' regDIF on.
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
#' @param estimator_history List to save EM iterations for supplemental EM algorithm.
#' @param estimator_limit Logical value indicating whether the EM algorithm reached
#' the maxit limit in the previous estimation round.
#' @param NA_cases Logical vector indicating if observation is missing.
#' @param exit_code Integer indicating if the model has converged properly.
#'
#' @return a \code{"list"} of matrices with unprocessed model estimates
#'
#' @keywords internal
#'
em_estimation <- function(p,
                          item_data,
                          pred_data,
                          prox_data,
                          mean_predictors,
                          var_predictors,
                          item_type,
                          theta,
                          pen_type,
                          tau_vec,
                          id_tau,
                          num_tau,
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
                          estimator_history,
                          estimator_limit,
                          NA_cases,
                          exit_code) {

  # Maximization and print settings.
  lastp <- p
  eps <- Inf
  iter <- 1
  models_to_fit <- ifelse(id_tau,num_tau,length(tau_vec))

  # Loop until convergence or maximum number of iterations.
  while(eps > final_control$tol && iter < final_control$maxit){


    # E-step: Evaluate Q function with current parameter estimates p.
    eout <- if(is.null(prox_data)) Estep(p,
                                         item_data,
                                         pred_data,
                                         item_type,
                                         mean_predictors,
                                         var_predictors,
                                         theta,
                                         samp_size,
                                         num_items,
                                         num_responses,
                                         adapt_quad,
                                         num_quad,
                                         get_eap = FALSE,
                                         NA_cases = NA_cases)



    if(optim_method == "MNR") {
      # M-step: Optimize parameters using multivariate NR.
      mout <- tryCatch(
        {
              Mstep_block(p,
                          item_data,
                          pred_data,
                          prox_data,
                          mean_predictors,
                          var_predictors,
                          eout,
                          item_type,
                          pen_type,
                          tau_vec[pen],
                          pen,
                          alpha,
                          gamma,
                          anchor,
                          final_control,
                          samp_size,
                          num_responses,
                          num_items,
                          num_quad,
                          num_predictors,
                          num_tau,
                          max_tau = FALSE)
        },
      error = function(e) {e; return(NULL)},
      warning = function(w) {} )

    } else if(optim_method == "UNR") {
      # M-step: Optimize parameters using one round of coordinate descent.
      mout <- tryCatch(
        {
          Mstep_cd(p,
                       item_data,
                       pred_data,
                       mean_predictors,
                       var_predictors,
                       eout,
                       item_type,
                       pen_type,
                       tau_vec[pen],
                       pen,
                       alpha,
                       gamma,
                       anchor,
                       final_control,
                       samp_size,
                       num_responses,
                       num_items,
                       num_quad,
                       num_predictors,
                       num_tau,
                       max_tau = FALSE)
    },
    error = function(e) {e; return(NULL)},
    warning = function(w) {})
    } else if(optim_method == "CD") {
      mout <- tryCatch(
        {
          Mstep_cd2(p,
                        item_data,
                        pred_data,
                        mean_predictors,
                        var_predictors,
                        eout,
                        item_type,
                        pen_type,
                        tau_vec[pen],
                        pen,
                        alpha,
                        gamma,
                        anchor,
                        final_control,
                        samp_size,
                        num_responses,
                        num_items,
                        num_quad,
                        num_predictors,
                        num_tau,
                        max_tau = FALSE)

        },
        error = function(e) {e; return(NULL)},
        warning = function(w) {})
    }


    if(is.null(mout)) {
      exit_code <- 4
      break
    }

    # Obtain parameter estimates.
    p <- mout$p

    # Update and check for convergence: Calculate the difference
    # in parameter estimates from current to previous.
    eps = sqrt(sum((unlist(p)-unlist(lastp))^2))

    # Save parameter estimates and observed log-likelihood for supplemental em.
    eout_obs_ll <- ifelse(is.null(eout), NA, eout$observed_ll)

    if(!is.null(eout)) {
      estimator_history[[pen]][,iter] <- c(unlist(p), eout_obs_ll)
    } else {
      estimator_history[[pen]][,iter] <- NA
    }


    # Add row for next EM step.
    if(eps > final_control$tol) {
      estimator_history[[pen]] <- cbind(estimator_history[[pen]],
                                 matrix(0,ncol=1,nrow=length(unlist(p))+1))
    }

    # Update parameter list.
    lastp <- p


    # Update the iteration number.
    iter = iter + 1
    if(iter == final_control$maxit) {
      warning("Iteration limit reached without convergence", call. = FALSE, immediate. = TRUE)
      estimator_limit <- T
      exit_code <- exit_code + 1
    }

      cat('\r', sprintf("Models Completed: %d of %d  Iteration: %d  Change: %f",
                        pen - 1,
                        models_to_fit,
                        iter,
                        round(eps, nchar(final_control$tol))))


    utils::flush.console()

    # Stop estimation if model would become under-identified because of tau
    # being too small.
    if(mout$under_identified) break
    # if(!is.null(prox_data)) break



  }

  if(exit_code == 4) return(NULL)

  # Get information criteria.
  infocrit <- information_criteria(eout,
                                   p,
                                   item_data,
                                   pred_data,
                                   prox_data,
                                   mean_predictors,
                                   var_predictors,
                                   item_type,
                                   gamma,
                                   samp_size,
                                   num_responses,
                                   num_items,
                                   num_quad)

  # Get EAP score.
  eout_eap <- if(is.null(prox_data)) {
    Estep(p,
          item_data,
          pred_data,
          item_type,
          mean_predictors,
          var_predictors,
          theta,
          samp_size,
          num_items,
          num_responses,
          adapt_quad,
          num_quad,
          get_eap = TRUE,
          NA_cases = NA_cases)
  } else {
    NULL
  }

  # Option to identify maximum value of tau which removes all DIF from model.
  if(id_tau) {

    # Final M-step.
    if(optim_method == "MNR") {
      max_tau <- Mstep_block(p,
                             item_data,
                             pred_data,
                             prox_data,
                             mean_predictors,
                             var_predictors,
                             eout,
                             item_type,
                             pen_type,
                             tau_vec[1],
                             pen,
                             alpha,
                             gamma,
                             anchor,
                             final_control,
                             samp_size,
                             num_responses,
                             num_items,
                             num_quad,
                             num_predictors,
                             num_tau,
                             max_tau = TRUE)
    } else if(optim_method == "UNR") {
      max_tau <- Mstep_cd(p,
                          item_data,
                          pred_data,
                          mean_predictors,
                          var_predictors,
                          eout,
                          item_type,
                          pen_type,
                          tau_vec[1],
                          pen,
                          alpha,
                          gamma,
                          anchor,
                          final_control,
                          samp_size,
                          num_responses,
                          num_items,
                          num_quad,
                          num_predictors,
                          num_tau,
                          max_tau = TRUE)
    } else if(optim_method == "CD") {
      max_tau <- Mstep_cd2(p,
                           item_data,
                           pred_data,
                           mean_predictors,
                           var_predictors,
                           eout,
                           item_type,
                           pen_type,
                           tau_vec[1],
                           pen,
                           alpha,
                           gamma,
                           anchor,
                           final_control,
                           samp_size,
                           num_responses,
                           num_items,
                           num_quad,
                           num_predictors,
                           num_tau,
                           max_tau = TRUE)
    }

  } else {
    max_tau <- NULL
  }

  # Return model results.
  return(list(p=p,
              complete_info=mout$inv_hess_diag,
              infocrit=infocrit,
              max_tau=max_tau,
              estimator_history=estimator_history,
              under_identified=mout$under_identified,
              estimator_limit=estimator_limit,
              eap=eout_eap,
              exit_code=exit_code))

}
