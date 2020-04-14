#######################
# Partial Derivatives #
#######################

#derivative for mean impact parameter
d_alpha <-
  function(p_impact,
           etable_all,
           theta,
           predictors,
           cov,
           samp_size,
           num_items,
           num_quadpts) {

  #get latent mean and variance vectors
  alpha <- predictors %*% p_impact[grep("g",names(p_impact),fixed=T)]
  phi <- exp(predictors %*% p_impact[grep("b",names(p_impact),fixed=T)])

  eta_d <- matrix(rep(predictors[,cov], num_quadpts), ncol = num_quadpts, nrow = samp_size)

  d1_trace <- t(sapply(1:samp_size, function(x){eta_d[x,]/phi[x]*(theta-alpha[x])}))
  d2_trace <- t(sapply(1:samp_size, function(x){-eta_d[x,]**2/phi[x]}))

  d1 <- sum(etable_all*d1_trace)
  d2 <- sum(etable_all*d2_trace)

  dlist <- list(d1,d2)

}

#derivative for variance impact parameter
d_phi <-
  function(p_impact,
           etable_all,
           theta,
           predictors,
           cov,
           samp_size,
           num_items,
           num_quadpts) {

  #get latent mean and variance vectors
  alpha <- predictors %*% p_impact[grep("g",names(p_impact),fixed=T)]
  phi <- exp(predictors %*% p_impact[grep("b",names(p_impact),fixed=T)])

  eta_d1 <- .5*sqrt(phi)*predictors[,cov]
  eta_d2 <- .25*sqrt(phi)*predictors[,cov]**2

  d1_trace <- t(sapply(1:samp_size, function(x) {eta_d1[x]*((theta-alpha[x])**2/phi[x]**(3/2) - 1/sqrt(phi[x]))}))
  d2_trace <- t(sapply(1:samp_size, function(x) {-2*eta_d2[x]*(phi[x]**(-3/2)*(theta-alpha[x])**2)}))

  d1 <- sum(etable_all*d1_trace)
  d2 <- sum(etable_all*d2_trace)

  dlist <- list(d1,d2)

}

#derivative for mu (mean) parameter
d_bernoulli <-
  function(parm,
           p_item,
           etable,
           theta,
           predictors,
           cov,
           samp_size,
           num_items,
           num_quadpts) {

  if(parm == "c0"){
    eta_d <- matrix(1, nrow = samp_size, ncol = num_quadpts)
  } else if(parm == "a0"){
    eta_d <- matrix(rep(matrix(theta, nrow = 1, ncol = length(theta)), samp_size), nrow = samp_size, ncol = length(theta), byrow = TRUE)
  } else if(parm == "c1"){
    eta_d <- matrix(rep(predictors[,cov], num_quadpts), ncol = num_quadpts, nrow = samp_size)
  } else if(parm == "a1"){
    eta_d <- matrix(rep(predictors[,cov], num_quadpts), ncol = num_quadpts, nrow = samp_size)*matrix(rep(matrix(theta, nrow = 1, ncol = length(theta)), samp_size), nrow = samp_size, ncol = length(theta), byrow = TRUE)
  }

  traceline <- bernoulli_traceline_pts(p_item,theta,predictors,samp_size,num_quadpts)

  d1 <- sum(traceline[[1]]*eta_d*etable[[2]]) + sum(-traceline[[2]]*eta_d*etable[[1]])
  d2 <- sum(-traceline[[2]]*traceline[[1]]*eta_d**2*etable[[1]]) + sum(-traceline[[2]]*traceline[[1]]*eta_d**2*etable[[2]])

  dlist <- list(d1,d2)

}

d_categorical <-
  function(parm,
           p_item,
           etable,
           theta,
           predictors,
           thr,
           cov,
           samp_size,
           num_responses_item,
           num_items,
           num_quadpts) {

  if(parm == "c0"){
    eta_d <- matrix(1, nrow = samp_size, ncol = num_quadpts)
  } else if(parm == "a0"){
    eta_d <- matrix(rep(matrix(theta, nrow = 1, ncol = length(theta)), samp_size), nrow = samp_size, ncol = length(theta), byrow = TRUE)
  } else if(parm == "c1"){
    eta_d <- matrix(rep(predictors[,cov], num_quadpts), ncol = num_quadpts, nrow = samp_size)
  } else if(parm == "a1"){
    eta_d <- matrix(rep(predictors[,cov], num_quadpts), ncol = num_quadpts, nrow = samp_size)*matrix(rep(matrix(theta, nrow = 1, ncol = length(theta)), samp_size), nrow = samp_size, ncol = length(theta), byrow = TRUE)
  }

  cum_traceline <- cumulative_traceline_pts(p_item,theta,predictors,samp_size,num_responses_item,num_quadpts)

  #non-threshold derivatives
  if(is.null(thr)){
    d1 <- eta_d*(-etable[[1]]*cum_traceline[[1]] +
                  etable[[num_responses_item]]*(1-cum_traceline[[num_responses_item-1]]))
    d2 <- eta_d**2*(-etable[[1]]*(cum_traceline[[1]]*(1-cum_traceline[[1]])) +
                     etable[[num_responses_item]]*(-cum_traceline[[num_responses_item-1]]*(1-cum_traceline[[num_responses_item-1]])))

    if(num_responses_item > 2) {
      for(i in 2:(num_responses_item-1)){
         #skip intermediate derivative calculations for constrained theshold on item 1
        d1 <- d1 + eta_d*etable[[i]]*((1-cum_traceline[[i]]) - cum_traceline[[i-1]])
        d2 <- d2 + eta_d**2*etable[[i]]*(cum_traceline[[i-1]]**2 + cum_traceline[[i]]**2 -
                                         cum_traceline[[i-1]] - cum_traceline[[i]])
      }
    }

    d1 <- sum(d1)
    d2 <- sum(d2)

    #threshold derivatives
  } else {
    cat_traceline <- categorical_traceline_pts(p_item,theta,predictors,samp_size,num_responses_item,num_quadpts)
    d1 <- sum(-etable[[thr]]*cum_traceline[[thr]]*(1-cum_traceline[[thr]])/cat_traceline[[thr]]) +
      sum(etable[[thr+1]]*cum_traceline[[thr]]*(1-cum_traceline[[thr]])/cat_traceline[[thr+1]])
    d2 <- sum(etable[[thr]]/cat_traceline[[thr]]*(cum_traceline[[thr]]*(1-cum_traceline[[thr]])**2 - cum_traceline[[thr]]**2*(1-cum_traceline[[thr]]) + cum_traceline[[thr]]**2*(1-cum_traceline[[thr]])**2/cat_traceline[[thr]])) -
      sum(etable[[thr+1]]/cat_traceline[[thr+1]]*(cum_traceline[[thr]]*(1-cum_traceline[[thr]])**2 - cum_traceline[[thr]]**2*(1-cum_traceline[[thr]]) - cum_traceline[[thr]]**2*(1-cum_traceline[[thr]])**2/cat_traceline[[thr+1]]))

  }

  dlist <- list(d1,d2)

}


d_mu_gaussian <-
  function(parm,
           p_item,
           etable,
           theta,
           responses_item,
           predictors,
           cov,
           samp_size,
           num_items,
           num_quadpts) {

  if(parm == "c0"){
    eta_d <- matrix(1, nrow = samp_size, ncol = num_quadpts)
  } else if(parm == "a0"){
    eta_d <- matrix(rep(matrix(theta, nrow = 1, ncol = length(theta)), samp_size), nrow = samp_size, ncol = length(theta), byrow = TRUE)
  } else if(parm == "c1"){
    eta_d <- matrix(rep(predictors[,cov], num_quadpts), ncol = num_quadpts, nrow = samp_size)
  } else if(parm == "a1"){
    eta_d <- matrix(rep(predictors[,cov], num_quadpts), ncol = num_quadpts, nrow = samp_size)*matrix(rep(matrix(theta, nrow = 1, ncol = length(theta)), samp_size), nrow = samp_size, ncol = length(theta), byrow = TRUE)
  }


  #get latent mean and variance vectors
  mu <- sapply(theta,function(x){(p_item[grep("c0",names(p_item),fixed=T)] + predictors %*% p_item[grep("c1",names(p_item),fixed=T)]) + (p_item[grep("a0",names(p_item),fixed=T)] + predictors %*% p_item[grep("a1",names(p_item),fixed=T)])*x})
  sigma <- sqrt(p_item[grep("s0",names(p_item))][1]*exp(predictors %*% p_item[grep("s1",names(p_item))]))


  d1_trace <- t(sapply(1:samp_size, function(x){eta_d[x,]/sigma[x]**2*(responses_item[x]-mu[x,])}))
  d2_trace <- t(sapply(1:samp_size, function(x){-eta_d[x,]**2/sigma[x]**2}))

  d1 <- sum(etable[[1]]*d1_trace)
  d2 <- sum(etable[[1]]*d2_trace)

  dlist <- list(d1,d2)

}

#derivative for sigma (variance) parameter


d_sigma_gaussian <-
  function(parm,
           p_item,
           etable,
           theta,
           responses_item,
           predictors,
           cov,
           samp_size,
           num_items,
           num_quadpts) {

  sigma <- sqrt(p_item[grep("s0",names(p_item))][1]*exp(predictors %*% p_item[grep("s1",names(p_item))]))
  mu <- sapply(theta,function(x){(p_item[grep("c0",names(p_item),fixed=T)] + predictors %*% p_item[grep("c1",names(p_item),fixed=T)]) + (p_item[grep("a0",names(p_item),fixed=T)] + predictors %*% p_item[grep("a1",names(p_item),fixed=T)])*x})

  if(parm == "s0"){
    eta_d1 <- sapply(1:samp_size, function(x) exp(predictors[x,] %*% p_item[grep("s1",names(p_item))])/(2*sigma[x]))
    eta_d2 <- sapply(1:samp_size, function(x) -exp(predictors[x,] %*% p_item[grep("s1",names(p_item))])**2/(4*sigma[x]**3))
  } else if(parm == "s1"){
    eta_d1 <- sapply(1:samp_size, function(x) sigma[x]*predictors[x,cov]/2)
    eta_d2 <- sapply(1:samp_size, function(x) sigma[x]*predictors[x,cov]**2/4)
  }


  d1_trace <- t(sapply(1:samp_size, function(x) {eta_d1[x]*((responses_item[x]-mu[x,])**2/sigma[x]**3 - 1/sigma[x])}))

  if(parm == "s0"){
    d2_trace <- t(sapply(1:samp_size, function(x) {eta_d1[x]**2*(1/sigma[x]**2 - 3*(responses_item[x]-mu[x,])**2/sigma[x]**4) +
                                                   eta_d2[x]*((responses_item[x]-mu[x,])**2/sigma[x]**3 - 1/sigma[x])}))
  } else if(parm == "s1"){
    d2_trace <- t(sapply(1:samp_size, function(x) {-2*eta_d2[x]*(sigma[x]**(-3)*(responses_item[x]-mu[x,])**2)}))
  }

  d1 <- sum(etable[[1]]*d1_trace)
  d2 <- sum(etable[[1]]*d2_trace)

  dlist <- list(d1,d2)

}



