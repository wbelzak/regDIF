###########################
# Log-likelihood - Impact #
###########################

ll.2pl.impact <-
  function(p_impact,nr,theta,covariates,samp_size,num_quadpts) {
    alpha <- covariates %*% p_impact[grep("g0",names(p_impact))]
    phi <- exp(covariates %*% p_impact[grep("b0",names(p_impact))])

    theta_scores2 <- t(sapply(1:samp_size, function(x){dnorm(theta, mean = alpha[x], sd = sqrt(phi[x]))}))

    ll_impact <- (-1)*(sum(nr*log(theta_scores)))
  }

