###########################
# Log-likelihood - Impact #
###########################

ll.2pl.impact <-
  function(p_impact,
           etable_all,
           theta,
           predictors,
           samp_size) {

    #get latent mean and variance vectors
    alpha <- predictors %*% p_impact[grep("g",names(p_impact),fixed=T)]
    phi <- exp(predictors %*% p_impact[grep("b",names(p_impact),fixed=T)])

    #get prior latent variable scores
    prior_scores <- t(sapply(1:samp_size, function(x){dnorm(theta, mean = alpha[x], sd = sqrt(phi[x]))/sum(dnorm(theta, mean = alpha[x], sd = sqrt(phi[x])))}))

    #log-likelihood for impact model
    ll_impact <- (-1)*(sum(etable_all*log(prior_scores)))
  }

