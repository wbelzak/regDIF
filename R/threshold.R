soft_threshold <- function(z,lambda){
  p_new <- sign(z)*max(abs(z) - lambda, 0)
  return(p_new)
}

firm_threshold <- function(z,lambda,gamma){
  if(abs(z) <= gamma*lambda){
    p_new <- (gamma/(gamma-1))*soft_threshold(z,lambda)
  }else{
    p_new <- z
  }
  return(p_new)
}
