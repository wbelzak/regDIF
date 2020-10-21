################
# Thresholding #
################

# Lasso penalty.
soft_threshold <-
  function(z,
           alpha,
           lambda) {

  p_new <- sign(z)*max(abs(z/(1+lambda*(1-alpha))) -
                         (lambda*alpha)/(1+lambda*(1-alpha)), 0)

  return(p_new)

}

# Mcp penalty.
firm_threshold <-
  function(z,
           alpha,
           lambda,
           gamma) {

  if(abs(z/(1+lambda*(1-alpha))) <= gamma*lambda){
    p_new <- (gamma/(gamma-1))*soft_threshold(z,alpha,lambda)
  }else{
    p_new <- z/(1+lambda*(1-alpha))
  }

  return(p_new)
}
