#' Expectation step.
#'
#' @param p List of parameters.
#' @param item.data Matrix or dataframe of item responses.
#' @param predictor.data Matrix or dataframe of DIF and/or impact predictors.
#' @param mean_predictors Possibly different matrix of predictors for the mean
#' impact equation.
#' @param var_predictors Possibly different matrix of predictors for the
#' variance impact equation.
#' @param samp_size Sample size in dataset.
#' @param num_items Number of items in dataset.
#' @param num_responses Number of responses for each item.
#' @param num_quadpts Number of quadrature points used for approximating the
#' latent variable.
#'
#' @NoRd
Estep <-
  function(p,
           item.data,
           predictor.data,
           mean_predictors,
           var_predictors,
           samp_size,
           num_items,
           num_responses,
           num_quadpts) {

  # Make space for the trace lines and the E-tables.
  itemtrace <- rep(list(NA),num_items)
  etable <- replicate(n=num_items,
                      matrix(0,nrow=samp_size,ncol=num_quadpts+1),
                      simplify = F)
  etable_all <- matrix(0,nrow=samp_size,ncol=num_quadpts)

  # Impact.
  alpha <- mean_predictors %*% p[[num_items+1]]
  phi <- exp(var_predictors %*% p[[num_items+2]])

  # Adaptive theta points.
  theta <- sapply(statmod::gauss.quad(n = num_quadpts, kind = "hermite")$nodes,
                  function(x) alpha + sqrt(2*phi)*x)

  # Compute the trace lines.
  for (item in 1:num_items) {
    if(num_responses[item] == 1) {
      itemtrace[[item]] <- gaussian_traceline_pts(p[[item]],
                                                  theta,
                                                  item.data[,item],
                                                  predictor.data,
                                                  samp_size,
                                                  num_quadpts)
    } else if (num_responses[item] == 2) {
      itemtrace[[item]] <- bernoulli_traceline_pts(p[[item]],
                                                   theta,
                                                   predictor.data,
                                                   samp_size,
                                                   num_quadpts)
    } else if (num_responses[item] > 2) {
      itemtrace[[item]] <- categorical_traceline_pts(p[[item]],
                                                     theta,
                                                     predictor.data,
                                                     samp_size,
                                                     num_responses[item],
                                                     num_quadpts)
    }
  }

  # Obtain E-tables.
  for(case in 1:samp_size) {

    # Adaptive qaudrature points.
    posterior <- dnorm(theta[case,], mean = alpha[case], sd = sqrt(phi[case]))

    # Within each response pattern, loop over items and compute posterior
    # probability of response pattern.
    for(item in 1:num_items) {
      x <- if(num_responses[item] == 1) {
             1
           } else {
             item.data[case,item]
           }
      if(!is.na(x)) posterior <- posterior*itemtrace[[item]][[x]][case,]
    }

    # Normalize posterior.
    marginal <- sum(posterior, na.rm = TRUE)
    if(marginal == 0) {marginal <- 1}
    posterior <- etable_all[case,] <- posterior/marginal

    # For individual i, add posterior to the E-tables depending on response.
    for(item in 1:num_items) {
      etable[[item]][case,] <- c(posterior, item.data[case,item])
    }

  }

  # List of E-tables to be used in Q function.
  elist <- list(etable = etable,
                etable_all = etable_all,
                theta = theta)


}

