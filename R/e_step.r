#####################################
# E-step: Evaluating the Q function #
#####################################

Estep_2pl <-
  function(p,
           responses,
           predictors,
           mean_predictors,
           var_predictors,
           theta,
           samp_size,
           num_items,
           num_responses,
           num_quadpts) { #p is parameters, responses is item responses, theta is values of latent variable

  #make space for the trace lines and the E-tables
  itemtrace <- rep(list(NA),num_items)
  etable <- replicate(n=num_items, matrix(0,nrow=samp_size,ncol=num_quadpts+1), simplify = F)
  etable_all <- matrix(0,nrow=samp_size,ncol=num_quadpts)

  #impact
  alpha <- mean_predictors %*% p[[num_items+1]]
  phi <- exp(var_predictors %*% p[[num_items+2]])

  #compute the trace lines
  for (item in 1:num_items) { #loop through items
    if(num_responses[item] == 1){
      itemtrace[[item]] <- gaussian_traceline_pts(p[[item]],theta,responses[,item],predictors,samp_size,num_quadpts)
    } else if (num_responses[item] == 2){
      itemtrace[[item]] <- bernoulli_traceline_pts(p[[item]],theta,predictors,samp_size,num_quadpts)
    } else if (num_responses[item] > 2){
      itemtrace[[item]] <- categorical_traceline_pts(p[[item]],theta,predictors,samp_size,num_responses[item],num_quadpts)
    }
  }

  #obtain E tables
  for(case in 1:samp_size) { #looping over samples

    #qaudrature points
    posterior <- dnorm(theta, mean = alpha[case], sd = sqrt(phi[case]))

    #within each response pattern, loop over items and compute posterior probability of response pattern
    for(item in 1:num_items) {
      x <- responses[case,item]
      if(!is.na(x)) posterior <- posterior*itemtrace[[item]][[x]][case,]
    }

    #normalize posterior to area equals number of persons with this response pattern
    marginal <- sum(posterior, na.rm = TRUE)
    if(marginal == 0) {marginal <- 1}
    posterior <- etable_all[case,] <- posterior/marginal

    #for individual i, add posterior to the r1 and r0 tables depending on response
    for(item in 1:num_items) { #within a person, loop over items
      x <- responses[case,item]
      etable[[item]][case,] <- c(posterior, x)
    }

  } #end loop over N

  #list of r tables to be used in Q function (ll.2pl function above)
  elist <- list(etable,etable_all)


} #end of E-step

