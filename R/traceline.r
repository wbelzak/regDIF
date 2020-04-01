#################
# 2PL Traceline #
#################
bernoulli_traceline_pts <-
  function(p_active,
           theta,
           predictors,
           samp_size,
           num_quadpts) {

    traceline <- replicate(n=2, matrix(0,nrow=samp_size,ncol=num_quadpts), simplify = F)

    traceline[[1]] <- sapply(theta,function(x){1-1/(1+exp(-((p_active[grep("c0",names(p_active),fixed=T)][1] + predictors %*% p_active[grep("c1",names(p_active),fixed=T)]) + (p_active[grep("a0",names(p_active),fixed=T)] + predictors %*% p_active[grep("a1",names(p_active),fixed=T)])*x)))})
    traceline[[2]] <- -traceline[[1]] + 1

    return(traceline)
  }

categorical_traceline_pts <-
  function(p_active,
           theta,
           predictors,
           samp_size,
           num_responses_item,
           num_quadpts) {

      #space for category traceline (y = c category)
      traceline <- replicate(n=num_responses_item, matrix(0,nrow=samp_size,ncol=num_quadpts), simplify = F)

      #for item response 1
      traceline[[1]] <- sapply(theta,function(x){1-1/(1+exp(-((p_active[grep("c0",names(p_active),fixed=T)][1] + predictors %*% p_active[grep("c1",names(p_active),fixed=T)]) + (p_active[grep("a0",names(p_active),fixed=T)] + predictors %*% p_active[grep("a1",names(p_active),fixed=T)])*x)))})
      #for item response 2
      if(num_responses_item > 2){
      traceline[[2]] <- sapply(theta,function(x){1/(1+exp(-((p_active[grep("c0",names(p_active),fixed=T)][1] + predictors %*% p_active[grep("c1",names(p_active),fixed=T)]) + (p_active[grep("a0",names(p_active),fixed=T)] + predictors %*% p_active[grep("a1",names(p_active),fixed=T)])*x)))}) -
                        sapply(theta,function(x){1/(1+exp(-((p_active[grep("c0",names(p_active),fixed=T)][1] - p_active[grep("c0",names(p_active),fixed=T)][2] + predictors %*% p_active[grep("c1",names(p_active),fixed=T)]) + (p_active[grep("a0",names(p_active),fixed=T)] + predictors %*% p_active[grep("a1",names(p_active),fixed=T)])*x)))})
        #for item responses 3 to J-1 (cycle through thresholds)
        if(num_responses_item > 3){
          for(thr in 3:(num_responses_item-1)){
            traceline[[thr]] <- sapply(theta,function(x){1/(1+exp(-((p_active[grep("c0",names(p_active),fixed=T)][1] - p_active[grep("c0",names(p_active),fixed=T)][thr-1] + predictors %*% p_active[grep("c1",names(p_active),fixed=T)]) + (p_active[grep("a0",names(p_active),fixed=T)] + predictors %*% p_active[grep("a1",names(p_active),fixed=T)])*x)))}) -
                                sapply(theta,function(x){1/(1+exp(-((p_active[grep("c0",names(p_active),fixed=T)][1] - p_active[grep("c0",names(p_active),fixed=T)][thr] + predictors %*% p_active[grep("c1",names(p_active),fixed=T)]) + (p_active[grep("a0",names(p_active),fixed=T)] + predictors %*% p_active[grep("a1",names(p_active),fixed=T)])*x)))})
          }
        }
      }
      #for item response J
      traceline[[num_responses_item]] <- sapply(theta,function(x){1/(1+exp(-((p_active[grep("c0",names(p_active),fixed=T)][1] - p_active[grep("c0",names(p_active),fixed=T)][num_responses_item-1] + predictors %*% p_active[grep("c1",names(p_active),fixed=T)]) + (p_active[grep("a0",names(p_active),fixed=T)] + predictors %*% p_active[grep("a1",names(p_active),fixed=T)])*x)))})

      return(traceline)
  }

cumulative_traceline_pts <-
  function(p_active,
           theta,
           predictors,
           samp_size,
           num_responses_item,
           num_quadpts) {

    #space for cumulative traceline (y >= c category)
    traceline <- replicate(n=(num_responses_item-1), matrix(0,nrow=samp_size,ncol=num_quadpts), simplify = F)

    #for item response 1
    traceline[[1]] <- sapply(theta,function(x){1/(1+exp(-((p_active[grep("c0",names(p_active),fixed=T)][1] + predictors %*% p_active[grep("c1",names(p_active),fixed=T)]) + (p_active[grep("a0",names(p_active),fixed=T)] + predictors %*% p_active[grep("a1",names(p_active),fixed=T)])*x)))})
    #for item response 2 to J
    if(num_responses_item > 2){
      for(thr in 2:(num_responses_item-1)){
        traceline[[thr]] <- sapply(theta,function(x){1/(1+exp(-((p_active[grep("c0",names(p_active),fixed=T)][1] - p_active[grep("c0",names(p_active),fixed=T)][thr] + predictors %*% p_active[grep("c1",names(p_active),fixed=T)]) + (p_active[grep("a0",names(p_active),fixed=T)] + predictors %*% p_active[grep("a1",names(p_active),fixed=T)])*x)))})
      }
    }

    return(traceline)
  }

gaussian_traceline_pts <-
  function(p_active,
           theta,
           responses_item,
           predictors,
           samp_size,
           num_quadpts) {

    # responses_item <- scale(responses_item)
    mu <- sapply(theta,function(x){(p_active[grep("c0",names(p_active),fixed=T)] + predictors %*% p_active[grep("c1",names(p_active),fixed=T)]) + (p_active[grep("a0",names(p_active),fixed=T)] + predictors %*% p_active[grep("a1",names(p_active),fixed=T)])*x})
    sigma <- sqrt(p_active[grep("s_",names(p_active))][1]*exp(predictors %*% p_active[grep("s(.*?)cov",names(p_active))]))

    traceline <- t(sapply(1:samp_size,function(x) dnorm(responses_item[x],mu[x,],sigma[x])))
    return(traceline)

  }



