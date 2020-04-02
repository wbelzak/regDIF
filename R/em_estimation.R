######################
# Reg-DIF Estimation #
######################

em_estimation <- function(p,
                          responses,
                          predictors,
                          theta,
                          itemtypes,
                          penalty,
                          lambda,
                          gamma,
                          pen,
                          anchor,
                          rasch,
                          final.control,
                          samp_size,
                          num_items,
                          num_responses,
                          num_predictors) {

  #Maximization settings
  lastp <- p
  eps <- Inf
  iter = 0

  #loop until convergence or maximum number of iterations
  while(eps > final.control$tol & iter < final.control$maxit){

    #E-step: Evaluate Q function with current parameter estimates p
    elist <- Estep_2pl(p,responses,predictors,theta,samp_size,num_items,num_responses,final.control$num_quadpts)

    #M-step: Optimize parameters
    p <- Mstep_2pl_dif(p,responses,predictors,elist,theta,itemtypes,penalty,lambda[pen],gamma,anchor,rasch,final.control$maxit,samp_size,num_responses,num_items,final.control$num_quadpts,num_predictors)

    #Update and check for convergence: Calculate the difference in parameter estimates from current to previous
    eps = sqrt(sum((unlist(p)-unlist(lastp))^2))

    #Update parameter list
    lastp <- p

    #update the iteration number
    iter = iter + 1
    if(iter == final.control$maxit) warning("EM iteration limit reached without convergence")

    cat('\r',"Models Completed:",pen-1,"of",length(lambda)," Iteration:",iter," Change:",round(eps,7)) #print information about optimization
    flush.console()

  } #End EM once converged or reached iteration limit

  return(list(elist,p))

}
