#######################
# Partial Derivatives #
#######################

#derivative for mu (mean) parameter
d_mu <-
  function(parm,
           p_item,
           etable,
           theta,
           responses_item,
           predictors,
           itemtypes,
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

    ##Categorical item responses##
    if(itemtypes == "categorical"){

      #non-threshold derivatives
      if(is.null(thr)){
        cat_traceline <- categorical_traceline_pts(p_item,theta,predictors,samp_size,num_responses_item,num_quadpts)
        cum_traceline <- cumulative_traceline_pts(p_item,theta,predictors,samp_size,num_responses_item,num_quadpts)
        d1 <- -eta_d*etable[[1]]*cum_traceline[[1]]*(1-cum_traceline[[1]])/cat_traceline[[1]] +
               eta_d*etable[[num_responses_item]]*cum_traceline[[num_responses_item-1]]*(1-cum_traceline[[num_responses_item-1]])/cat_traceline[[num_responses_item]]
        d2 <- -eta_d**2*etable[[1]]/cat_traceline[[1]]**2*((cum_traceline[[1]]*(1-cum_traceline[[1]])**2 - cum_traceline[[1]]**2*(1-cum_traceline[[1]]))*cat_traceline[[1]] + cum_traceline[[1]]**2*(1-cum_traceline[[1]])**2) +
               eta_d**2*etable[[num_responses_item]]/cat_traceline[[num_responses_item]]**2*((cum_traceline[[num_responses_item-1]]*(1-cum_traceline[[num_responses_item-1]])**2 - cum_traceline[[num_responses_item-1]]**2*(1-cum_traceline[[num_responses_item-1]]))*cat_traceline[[num_responses_item]] - cum_traceline[[num_responses_item-1]]**2*(1-cum_traceline[[num_responses_item-1]])**2)
        for(i in 2:(num_responses_item-1)){
          if(all(cat_traceline[[i]] == 0) | num_responses_item == 2) {next} #skip intermediate derivative calculations for constrained theshold on item 1
          d1 <- d1 + eta_d*etable[[i]]/cat_traceline[[i]]*(cum_traceline[[i-1]]*(1-cum_traceline[[i-1]]) - cum_traceline[[i]]*(1-cum_traceline[[i]]))
          d2 <- d2 + eta_d**2*etable[[i]]/cat_traceline[[i]]**2*((cum_traceline[[i-1]]*(1-cum_traceline[[i-1]])**2 - cum_traceline[[i-1]]**2*(1-cum_traceline[[i-1]]) -
                                                                  cum_traceline[[i]]*(1-cum_traceline[[i]])**2 + cum_traceline[[i]]**2*(1-cum_traceline[[i]]))*cat_traceline[[i]] -
                                                                 (cum_traceline[[i-1]]*(1-cum_traceline[[i-1]]) - cum_traceline[[i]]*(1-cum_traceline[[i]]))**2)
        }
        d1 <- sum(d1)
        d2 <- sum(d2)

      #threshold derivatives
      } else {
        cat_traceline <- categorical_traceline_pts(p_item,theta,predictors,samp_size,num_responses_item,num_quadpts)
        cum_traceline <- cumulative_traceline_pts(p_item,theta,predictors,samp_size,num_responses_item,num_quadpts)
        d1 <- sum(-etable[[thr]]*cum_traceline[[thr]]*(1-cum_traceline[[thr]])/cat_traceline[[thr]]) +
              sum(etable[[thr+1]]*cum_traceline[[thr]]*(1-cum_traceline[[thr]])/cat_traceline[[thr+1]])
        d2 <- sum(etable[[thr]]/cat_traceline[[thr]]*(cum_traceline[[thr]]*(1-cum_traceline[[thr]])**2 - cum_traceline[[thr]]**2*(1-cum_traceline[[thr]]) + cum_traceline[[thr]]**2*(1-cum_traceline[[thr]])**2/cat_traceline[[thr]])) -
              sum(etable[[thr+1]]/cat_traceline[[thr+1]]*(cum_traceline[[thr]]*(1-cum_traceline[[thr]])**2 - cum_traceline[[thr]]**2*(1-cum_traceline[[thr]]) - cum_traceline[[thr]]**2*(1-cum_traceline[[thr]])**2/cat_traceline[[thr+1]]))
      }

   ##Continuous item responses##
   } else if(itemtypes == "continuous"){

     # responses_item <- scale(responses_item)
     mu <- sapply(theta,function(x){(p_item[grep("c0",names(p_item),fixed=T)] + predictors %*% p_item[grep("c1",names(p_item),fixed=T)]) + (p_item[grep("a0",names(p_item),fixed=T)] + predictors %*% p_item[grep("a1",names(p_item),fixed=T)])*x})
     sigma <- sqrt(p_item[grep("s_",names(p_item))][1]*exp(predictors %*% p_item[grep("s(.*?)cov",names(p_item))]))
     d1_trace <- t(sapply(1:samp_size, function(x) eta_d[x,]*(responses_item[x]-mu[x,])/sigma[x]**2))
     d2_trace <- t(sapply(1:samp_size, function(x) -eta_d[x,]**2/sigma[x]**2))

     d1 <- sum(etable[[1]]*d1_trace)
     d2 <- sum(etable[[1]]*d2_trace)

    }
    dlist <- list(d1,d2)
  }

#derivative for sigma (variance) parameter
d_sigma <-
  function(parm,
           p_item,
           etable,
           theta,
           responses_item,
           predictors,
           itemtypes,
           cov,
           samp_size,
           num_items,
           num_quadpts) {

    sigma <- sqrt(p_item[grep("s_",names(p_item))][1]*exp(predictors %*% p_item[grep("s(.*?)cov",names(p_item))]))
    # responses_item <- scale(responses_item)
    mu <- sapply(theta,function(x){(p_item[grep("c0",names(p_item),fixed=T)] + predictors %*% p_item[grep("c1",names(p_item),fixed=T)]) + (p_item[grep("a0",names(p_item),fixed=T)] + predictors %*% p_item[grep("a1",names(p_item),fixed=T)])*x})

    if(parm == "s0"){
      eta_d1 <- sapply(1:samp_size, function(x) (sigma[x]**2/as.numeric(p_item[grep("s_",names(p_item))][1]))/(2*as.numeric(p_item[grep("s_",names(p_item))][1])/sigma[x]))
      eta_d2 <- sapply(1:samp_size, function(x) -sigma[x]**2/as.numeric(p_item[grep("s_",names(p_item))][1])/4*as.numeric(p_item[grep("s_",names(p_item))][1])**(3/2))
    } else if(parm == "s1"){
      eta_d1 <- sapply(1:samp_size, function(x) sigma[x]*predictors[x,cov]/2)
      eta_d2 <- sapply(1:samp_size, function(x) sigma[x]*predictors[x,cov]**2/4)
    }

    d1_trace <- t(sapply(1:samp_size, function(x) eta_d1[x]*((responses_item[x]-mu[x,])**2/sigma[x]**3 - 1/sigma[x])))
    d2_trace <- t(sapply(1:samp_size, function(x) eta_d1[x]**2*((1/sigma[x]**2) - 3*(responses_item[x]-mu[x,])**2/sigma[x]**4) + eta_d2[x]*((responses_item[x]-mu[x,])**2/sigma[x]**3 - 1/sigma[x])))

    d1 <- sum(etable[[1]]*d1_trace)
    d2 <- sum(etable[[1]]*d2_trace)

    dlist <- list(d1,d2)

  }



