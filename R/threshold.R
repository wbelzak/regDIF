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

# Lasso group penalty.
grp_soft_threshold <-
  function(z,
           tau) {

    l2_norm_z <- sqrt(sum(z**2))
    p_new <-
      if(l2_norm_z > tau) {
        (l2_norm_z - tau)*(z/l2_norm_z)
      } else if (l2_norm_z <= tau) {
        c(0,0)
      }

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

# Mcp penalty.
grp_firm_threshold <-
  function(z,
           tau,
           gamma) {

    l2_norm_z <- sqrt(sum(z**2))
    if(l2_norm_z <= gamma*tau){
      p_new <- (gamma/(gamma-1))*grp_soft_threshold(z,tau)
    }else{
      p_new <- z
    }

    return(p_new)
  }
