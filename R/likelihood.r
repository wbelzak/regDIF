###########################
# Log-likelihood - Impact #
###########################

ll.2pl.impact <-
  function(p_impact,nr,theta,covariates,samp_size,num_quadpts) {
    alpha <- covariates %*% p_impact[grep("g0",names(p_impact))]
    phi <- exp(covariates %*% p_impact[grep("b0",names(p_impact))])

    theta_scores <- matrix(0,nrow=samp_size,ncol=num_quadpts)
    for(case in 1:samp_size){
      theta_scores[case,] <- dnorm(theta, mean = alpha[case], sd = sqrt(phi[case]))
    }

    ll_impact <- (-1)*(sum(nr*log(theta_scores)))
  }

