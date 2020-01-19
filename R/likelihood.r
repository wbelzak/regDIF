###########################
# Log-likelihood - Impact #
###########################

ll.2pl.impact <-
  function(p_impact,nr,theta,covariates,samp_size,num_quadpts) {
    alpha <- covariates %*% p_impact[grep("g0",names(p_impact),fixed=T)]
    phi <- exp(covariates %*% p_impact[grep("b0",names(p_impact),fixed=T)])

    theta_scores <- t(sapply(1:samp_size, function(x){dnorm(theta, mean = alpha[x], sd = sqrt(phi[x]))}))

    ll_impact <- (-1)*(sum(nr*log(theta_scores)))
  }

