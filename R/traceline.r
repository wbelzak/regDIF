##################
# 2PL Tracelines #
##################

bernoulli_traceline_pts <-
  function(p_active,
           theta,
           predictor.data,
           alpha,
           phi,
           samp_size,
           num_quadpts) {

    c0_parms <- grepl("c0",names(p_active),fixed=T)
    c1_parms <- grepl("c1",names(p_active),fixed=T)
    a0_parms <- grepl("a0",names(p_active),fixed=T)
    a1_parms <- grepl("a1",names(p_active),fixed=T)

    traceline0 <- apply(theta, 2, function(x){1-1/(1+exp(-((p_active[c0_parms] + predictor.data %*% p_active[c1_parms]) + (p_active[a0_parms] + predictor.data %*% p_active[a1_parms])*x)))})
    traceline1 <- -traceline0 + 1

    return(list(traceline0,traceline1))

  }

categorical_traceline_pts <-
  function(p_active,
           theta,
           predictor.data,
           samp_size,
           num_responses_item,
           num_quadpts) {

  #space for category traceline (y = c category)
  traceline <- replicate(n=num_responses_item, matrix(0,nrow=samp_size,ncol=num_quadpts), simplify = F)


  #for item response 1
  traceline[[1]] <- sapply(theta,function(x){1-1/(1+exp(-((p_active[grep("c0",names(p_active),fixed=T)][1] + predictor.data %*% p_active[grep("c1",names(p_active),fixed=T)]) + (p_active[grep("a0",names(p_active),fixed=T)] + predictor.data %*% p_active[grep("a1",names(p_active),fixed=T)])*x)))})
  #for item response 2
  if(num_responses_item > 2){
  traceline[[2]] <- sapply(theta,function(x){1/(1+exp(-((p_active[grep("c0",names(p_active),fixed=T)][1] + predictor.data %*% p_active[grep("c1",names(p_active),fixed=T)]) + (p_active[grep("a0",names(p_active),fixed=T)] + predictor.data %*% p_active[grep("a1",names(p_active),fixed=T)])*x)))}) -
                    sapply(theta,function(x){1/(1+exp(-((p_active[grep("c0",names(p_active),fixed=T)][1] - p_active[grep("c0",names(p_active),fixed=T)][2] + predictor.data %*% p_active[grep("c1",names(p_active),fixed=T)]) + (p_active[grep("a0",names(p_active),fixed=T)] + predictor.data %*% p_active[grep("a1",names(p_active),fixed=T)])*x)))})
    #for item responses 3 to J-1 (cycle through thresholds)
    if(num_responses_item > 3){
      for(thr in 3:(num_responses_item-1)){
        traceline[[thr]] <- sapply(theta,function(x){1/(1+exp(-((p_active[grep("c0",names(p_active),fixed=T)][1] - p_active[grep("c0",names(p_active),fixed=T)][thr-1] + predictor.data %*% p_active[grep("c1",names(p_active),fixed=T)]) + (p_active[grep("a0",names(p_active),fixed=T)] + predictor.data %*% p_active[grep("a1",names(p_active),fixed=T)])*x)))}) -
                            sapply(theta,function(x){1/(1+exp(-((p_active[grep("c0",names(p_active),fixed=T)][1] - p_active[grep("c0",names(p_active),fixed=T)][thr] + predictor.data %*% p_active[grep("c1",names(p_active),fixed=T)]) + (p_active[grep("a0",names(p_active),fixed=T)] + predictor.data %*% p_active[grep("a1",names(p_active),fixed=T)])*x)))})
      }
    }
  }
  #for item response J
  traceline[[num_responses_item]] <- sapply(theta,function(x){1/(1+exp(-((p_active[grep("c0",names(p_active),fixed=T)][1] - p_active[grep("c0",names(p_active),fixed=T)][num_responses_item-1] + predictor.data %*% p_active[grep("c1",names(p_active),fixed=T)]) + (p_active[grep("a0",names(p_active),fixed=T)] + predictor.data %*% p_active[grep("a1",names(p_active),fixed=T)])*x)))})

  return(traceline)

}

cumulative_traceline_pts <-
  function(p_active,
           theta,
           predictor.data,
           samp_size,
           num_responses_item,
           num_quadpts) {

  #space for cumulative traceline (y >= c category)
  traceline <- replicate(n=(num_responses_item-1), matrix(0,nrow=samp_size,ncol=num_quadpts), simplify = F)

  #for item response 1
  traceline[[1]] <- sapply(theta,function(x){1/(1+exp(-((p_active[grep("c0",names(p_active),fixed=T)][1] + predictor.data %*% p_active[grep("c1",names(p_active),fixed=T)]) + (p_active[grep("a0",names(p_active),fixed=T)] + predictor.data %*% p_active[grep("a1",names(p_active),fixed=T)])*x)))})
  #for item response 2 to J
  if(num_responses_item > 2){
    for(thr in 2:(num_responses_item-1)){
      traceline[[thr]] <- sapply(theta,function(x){1/(1+exp(-((p_active[grep("c0",names(p_active),fixed=T)][1] - p_active[grep("c0",names(p_active),fixed=T)][thr] + predictor.data %*% p_active[grep("c1",names(p_active),fixed=T)]) + (p_active[grep("a0",names(p_active),fixed=T)] + predictor.data %*% p_active[grep("a1",names(p_active),fixed=T)])*x)))})
    }
  }

  return(traceline)

}

gaussian_traceline_pts <-
  function(p_active,
           theta,
           responses_item,
           predictor.data,
           samp_size,
           num_quadpts) {

  # responses_item <- scale(responses_item)
  mu <- sapply(theta,function(x){(p_active[grep("c0",names(p_active),fixed=T)] + predictor.data %*% p_active[grep("c1",names(p_active),fixed=T)]) + (p_active[grep("a0",names(p_active),fixed=T)] + predictor.data %*% p_active[grep("a1",names(p_active),fixed=T)])*x})
  sigma <- sqrt(p_active[grep("s0",names(p_active))][1]*exp(predictor.data %*% p_active[grep("s1",names(p_active))]))

  traceline <- t(sapply(1:samp_size,function(x) dnorm(responses_item[x],mu[x,],sigma[x])))

  return(list(traceline))

}



