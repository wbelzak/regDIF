#####################################
# E-step: Evaluating the Q function #
#####################################

Estep_2pl <-
  function(p,
           item.data,
           predictor.data,
           mean_predictors,
           var_predictors,
           samp_size,
           num_items,
           num_responses,
           num_quadpts) {

  #make space for the trace lines and the E-tables
  itemtrace <- rep(list(NA),num_items)
  etable <- replicate(n=num_items, matrix(0,nrow=samp_size,ncol=num_quadpts+1), simplify = F)
  etable_all <- matrix(0,nrow=samp_size,ncol=num_quadpts)

  #impact
  alpha <- mean_predictors %*% p[[num_items+1]]
  phi <- exp(var_predictors %*% p[[num_items+2]])

  #adaptive theta
  theta <- sapply(statmod::gauss.quad(n = num_quadpts, kind = "hermite")$nodes,
                  function(x) alpha + sqrt(2*phi)*x)

  #compute the trace lines
  for (item in 1:num_items) { #loop through items
    if(num_responses[item] == 1){
      itemtrace[[item]] <- gaussian_traceline_pts(p[[item]],theta,item.data[,item],predictor.data,samp_size,num_quadpts)
    } else if (num_responses[item] == 2){
      itemtrace[[item]] <- bernoulli_traceline_pts(p[[item]],theta,predictor.data,samp_size,num_quadpts)
    } else if (num_responses[item] > 2){
      itemtrace[[item]] <- categorical_traceline_pts(p[[item]],theta,predictor.data,samp_size,num_responses[item],num_quadpts)
    }
  }

  #obtain E tables
  for(case in 1:samp_size) { #looping over samples

    #qaudrature points
    posterior <- dnorm(theta[case,], mean = alpha[case], sd = sqrt(phi[case]))

    #within each response pattern, loop over items and compute posterior probability of response pattern
    for(item in 1:num_items) {
      x <- if(num_responses[item] == 1) {1} else {item.data[case,item]}
      if(!is.na(x)) posterior <- posterior*itemtrace[[item]][[x]][case,]
    }

    #normalize posterior
    marginal <- sum(posterior, na.rm = TRUE)
    if(marginal == 0) {marginal <- 1}
    posterior <- etable_all[case,] <- posterior/marginal

    #for individual i, add posterior to the e-tables depending on response
    for(item in 1:num_items) { #within a person, loop over items
      x <- item.data[case,item]
      etable[[item]][case,] <- c(posterior, x)
    }

  } #end loop over N

  #list of r tables to be used in Q function (ll.2pl function above)
  elist <- list(etable = etable,
                etable_all = etable_all,
                theta = theta)


} #end of E-step

