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


########################
# Log-likelihood - DIF #
########################

ll.2pl.dif <-
  function(p_active,r1,r0,theta,covariates,tau,samp_size,num_quadpts) { #p2 is current estimate of parameters, r1 is the conditional expected proportion of individuals that endorse item i, r0 is the conditional expected proportion of individuals that do not endorse item i

    itemtrace <- matrix(0,nrow=samp_size,ncol=num_quadpts)
    for(i in 1:num_quadpts){
      itemtrace[,i] <- trace.line.pts(p_active,theta[i],covariates)
    }

    ll_dif <- (-1)*(sum(r1*(log(itemtrace)), na.rm = TRUE) + sum(r0*(log(1.0-itemtrace)), na.rm = TRUE))
    pen <- (-1)*(tau*sum(c(abs(p_active[grep("c1_itm",names(p_active))]), abs(p_active[grep("a1_itm",names(p_active))])), na.rm = TRUE))
    ll_dif_pen <- ll_dif - pen #Q function we want to minimize
  }
