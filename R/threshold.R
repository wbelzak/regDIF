################
# Thresholding #
################

# Lasso penalty.
soft_threshold <-
  function(z,
           alpha,
           tau) {

  p_new <- sign(z)*max(abs(z/(1+tau*(1-alpha))) -
                         (tau*alpha)/(1+tau*(1-alpha)), 0)

  return(p_new)

}

# Mcp penalty.
firm_threshold <-
  function(z,
           alpha,
           tau,
           gamma) {

  if(abs(z/(1+tau*(1-alpha))) <= gamma*tau){
    p_new <- (gamma/(gamma-1))*soft_threshold(z,alpha,tau)
  }else{
    p_new <- z/(1+tau*(1-alpha))
  }

  return(p_new)
}
