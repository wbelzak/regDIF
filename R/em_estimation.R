######################
# Reg-DIF Estimation #
######################

em_estimation <- function(p,
                          responses,
                          predictors,
                          mean_predictors,
                          var_predictors,
                          theta,
                          itemtypes,
                          penalty,
                          lambda,
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
                          num_quadpts) {

   # p <- data_scrub$p; responses <- data_scrub$responses;predictors <- data_scrub$predictors;mean_predictors <- data_scrub$mean_predictors;var_predictors <- data_scrub$var_predictors;theta <- data_scrub$theta;itemtypes <- data_scrub$itemtypes;lambda <- data_scrub$lambda; final.control <- data_scrub$final.control;samp_size <- data_scrub$samp_size;num_items <- data_scrub$num_items;num_responses <- data_scrub$num_responses;num_predictors <- data_scrub$num_predictors;num_quadpts <- data_scrub$num_quadpts
  #Maximization settings
  lastp <- p
  eps <- Inf
  iter = 0

  #loop until convergence or maximum number of iterations
  while(eps > final.control$tol & iter < final.control$maxit){

    #E-step: Evaluate Q function with current parameter estimates p
    elist <- Estep_2pl(p,responses,predictors,mean_predictors,var_predictors,theta,samp_size,num_items,num_responses,num_quadpts)

    #M-step: Optimize parameters
    p <- Mstep_2pl_dif(p,responses,predictors,mean_predictors,var_predictors,elist,theta,itemtypes,penalty,lambda[pen],alpha,gamma,anchor,rasch,samp_size,num_responses,num_items,num_quadpts,num_predictors)

    # p <- em_step(p,theta,responses,predictors,itemtypes,penalty,lambda,pen,alpha,gamma,anchor,rasch,samp_size,num_items,num_responses,num_quadpts,num_predictors)

    #Update and check for convergence: Calculate the difference in parameter estimates from current to previous
    eps = sqrt(sum((unlist(p)-unlist(lastp))^2))

    #Update parameter list
    lastp <- p

    #update the iteration number
    iter = iter + 1
    if(iter == final.control$maxit) warning("EM iteration limit reached without convergence")

    cat('\r',sprintf("Models Completed: %d of %d  Iteration: %d  Change: %f", pen-1, length(lambda), iter, round(eps,6))) #print information about optimization

    utils::flush.console()

  } #End EM once converged or reached iteration limit

  #get information criteria
  infocrit <- information_criteria(elist,p,responses,predictors,mean_predictors,var_predictors,theta,lambda[pen],samp_size,num_responses,num_items,num_quadpts)

  return(list(p,infocrit))

}
