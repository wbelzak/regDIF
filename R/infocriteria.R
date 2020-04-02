########################
# Information Criteria #
########################

information_criteria <-
  function(elist,
           p,
           theta,
           predictors,
           lambda,
           pen,
           samp_size,
           num_responses,
           num_items,
           num_quadpts) {

  ll_dif <- 0
  for (item in 1:num_items) { #loop through items

    #obtain etables for each response category
    etable <- replicate(n=num_responses[item], elist[[1]][[item]], simplify = F)
    for(resp in 1:num_responses[item]){
      etable[[resp]][which(!(etable[[resp]][,ncol(etable[[resp]])] == resp)),] <- 0
    }
    etable <- lapply(etable, function(x) x[,1:num_quadpts])

    #get item parameters
    p_item <- p[[item]]

    #compute negative log-likelihood values
    if(num_responses[item] == 1){
      itemtrace <- gaussian_traceline_pts(p[[item]],theta,responses[,item],predictors,samp_size,num_quadpts)
      ll_dif_item <- -1*sum(etable[[1]]*log(itemtrace), na.rm = TRUE)
    } else if (num_responses[item] == 2){
      itemtrace <- bernoulli_traceline_pts(p[[item]],theta,predictors,samp_size,num_quadpts)
      ll_dif_item <- -1*(sum(etable[[1]]*log(itemtrace[[1]]), na.rm = TRUE) + sum(etable[[2]]*log(itemtrace[[2]]), na.rm = TRUE))
    } else if (num_responses[item] > 2){
      itemtrace <- categorical_traceline_pts(p[[item]],theta,predictors,samp_size,num_responses[item],num_quadpts)
      ll_dif_item <- 0
      for(resp in 1:num_responses[item]){
      ll_dif_item <- ll_dif_item + -1*sum(etable[[resp]]*ifelse(all(itemtrace[[item]][[resp]] == 0),0,itemtrace[[resp]]), na.rm = TRUE)
      }
    }

    #subtract out penalty and add negative log-likelihood over all items
    current_pen <- (lambda[pen]*sum(c(abs(p_item[grep("c1_itm",names(p_item))]), abs(p_item[grep("a1_itm",names(p_item))])), na.rm = TRUE))
    ll_dif <- ll_dif + ll_dif_item - current_pen #Q function we want to minimize
  }

  #obtain likelihood value for latent variable model
  p_impact <- c(p[[num_items+1]],p[[num_items+2]])
  alpha <- predictors %*% p_impact[grep("g",names(p_impact),fixed=T)]
  phi <- exp(predictors %*% p_impact[grep("b",names(p_impact),fixed=T)])
  prior_scores <- t(sapply(1:samp_size, function(x){dnorm(theta, mean = alpha[x], sd = sqrt(phi[x]))}))
  ll_impact <- -1*sum(elist[[2]]*log(prior_scores))

  #remove dif parameters that equal zero from information criterion calculation
  p2 <- unlist(p)
  p2_base <- p2[c(grep("c0",names(p2)),grep("a0",names(p2)))]
  p2_cov <- p2[c(grep("c1",names(p2)),grep("a1",names(p2)))]
  p2_cov <- p2_cov[p2_cov != 0]
  p2 <- c(p2_base,p2_cov,p_impact)
  ll <- ll_impact + ll_dif

  #compute aic and bic
  aic <- 2*length(p2) + 2*ll
  bic <- log(samp_size)*length(p2) + 2*ll

  return(c(aic,bic))

}
