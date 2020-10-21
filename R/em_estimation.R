######################
# Reg-DIF Estimation #
######################

em_estimation <- function(p,
                          item.data,
                          predictor.data,
                          mean_predictors,
                          var_predictors,
                          item.type,
                          penalty,
                          tau,
                          alpha,
                          gamma,
                          pen,
                          anchor,
                          rasch,
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
    elist <- Estep_2pl(p,
                       item.data,
                       predictor.data,
                       mean_predictors,
                       var_predictors,
                       samp_size,
                       num_items,
                       num_responses,
                       num_quadpts)

    # M-step: Optimize parameters.
    p <- Mstep_2pl_dif(p,
                       item.data,
                       predictor.data,
                       mean_predictors,
                       var_predictors,
                       elist,
                       item.type,
                       penalty,
                       tau[pen],
                       alpha,
                       gamma,
                       anchor,
                       rasch,
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
                                   penalty,
                                   samp_size,
                                   num_responses,
                                   num_items,
                                   num_quadpts)

  return(list(p,infocrit))

}
