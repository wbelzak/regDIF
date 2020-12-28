#' Expectation step.
#'
#' @param p List of parameters.
#' @param item.data Matrix or dataframe of item responses.
#' @param pred.data Matrix or dataframe of DIF and/or impact predictors.
#' @param mean_predictors Possibly different matrix of predictors for the mean
#' impact equation.
#' @param var_predictors Possibly different matrix of predictors for the
#' variance impact equation.
#' @param samp_size Sample size in dataset.
#' @param num_items Number of items in dataset.
#' @param num_responses Number of responses for each item.
#' @param num.quad Number of quadrature points used for approximating the
#' latent variable.
#'
#' @keywords internal
#'
Estep <-
  function(p,
           item.data,
           pred.data,
           mean_predictors,
           var_predictors,
           samp_size,
           num_items,
           num_responses,
           adapt.quad,
           num.quad) {

  # Make space for the trace lines and the E-tables.
  itemtrace <- rep(list(NA),num_items)
  etable <- replicate(n=num_items,
                      matrix(0,nrow=samp_size,ncol=num.quad+1),
                      simplify = F)
  etable_all <- matrix(0,nrow=samp_size,ncol=num.quad)

  # Impact.
  alpha <- mean_predictors %*% p[[num_items+1]]
  phi <- exp(var_predictors %*% p[[num_items+2]])

  # Adaptive theta points.
  if(adapt.quad == TRUE) {
    theta <- sapply(statmod::gauss.quad(n = num.quad,
                                        kind = "hermite")$nodes,
                    function(x) alpha + sqrt(2*phi)*x)
  } else {
    theta <- matrix(seq(-6, 6, length.out = num.quad), nrow = samp_size,
                    ncol = num.quad, byrow = T)
  }


  # Compute the trace lines.
  for (item in 1:num_items) {
    if(num_responses[item] == 1) {
      itemtrace[[item]] <- gaussian_traceline_pts(p[[item]],
                                                  theta,
                                                  item.data[,item],
                                                  pred.data,
                                                  samp_size,
                                                  num.quad)
    } else if (num_responses[item] == 2) {
      itemtrace[[item]] <- bernoulli_traceline_cpp(p[[item]],
                                                   theta,
                                                   pred.data,
                                                   samp_size,
                                                   num.quad)
    } else if (num_responses[item] > 2) {
      itemtrace[[item]] <- categorical_traceline_cpp(p[[item]],
                                                     theta,
                                                     pred.data,
                                                     samp_size,
                                                     num_responses[item],
                                                     num.quad)
    }
  }

  # Obtain E-tables.
  for(case in 1:samp_size) {

    # Get weights for computing adaptive or fixed-point quadrature.
    posterior <- dnorm(theta[case,],
                       mean = alpha[case],
                       sd = sqrt(phi[case]))

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

