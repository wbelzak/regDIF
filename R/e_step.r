#' Expectation step.
#'
#' @param p List of parameters.
#' @param item_data Matrix or dataframe of item responses.
#' @param pred_data Matrix or dataframe of DIF and/or impact predictors.
#' @param item_type Vector of character values indicating the item type.
#' @param mean_predictors Possibly different matrix of predictors for the mean
#' impact equation.
#' @param var_predictors Possibly different matrix of predictors for the
#' variance impact equation.
#' @param theta Vector of fixed quadrature points.
#' @param samp_size Sample size in dataset.
#' @param num_items Number of items in dataset.
#' @param num_responses Number of responses for each item.
#' @param num_quad Number of quadrature points used for approximating the
#' latent variable.
#' @param get_eap Logical indicating whether to compute EAP scores.
#' @param NA_cases Logical vector indicating missing observations.
#'
#' @return a \code{"list"} of posterior values from the expectation step
#'
#' @keywords internal
#'
Estep <-
  function(p,
           item_data,
           pred_data,
           item_type,
           mean_predictors,
           var_predictors,
           theta,
           samp_size,
           num_items,
           num_responses,
           adapt_quad,
           num_quad,
           get_eap,
           NA_cases) {

    # Make space for the trace lines and the E-tables.
    observed_ll <- 0
    itemtrace <- rep(list(NA),num_items)
    etable <- matrix(0, nrow = samp_size, ncol = num_quad)
    eap_scores <- if(get_eap) matrix(NA, nrow = length(NA_cases), ncol = 1)
    eap_sd <- if(get_eap) matrix(NA, nrow = length(NA_cases), ncol = 1)

    # Impact.
    alpha <- mean_predictors %*% p[[num_items+1]]
    phi <- exp(var_predictors %*% p[[num_items+2]])

    if(adapt_quad == TRUE) {
      theta <- mean(alpha) +
        sqrt(2*mean(phi))*statmod::gauss.quad(n = num_quad,
                                              kind = "hermite")$nodes
    }

    # Compute the trace lines.
    for (item in 1:num_items) {
      if(item_type[item] == "cfa") {
        itemtrace[[item]] <- gaussian_traceline_pts(p[[item]],
                                                    theta,
                                                    item_data[,item],
                                                    pred_data,
                                                    samp_size)
      } else if (num_responses[item] == 2) {
        itemtrace[[item]] <- bernoulli_traceline_pts(p[[item]],
                                                     theta,
                                                     pred_data,
                                                     samp_size)
      } else {
        itemtrace[[item]] <- cumulative_traceline_pts(p[[item]],
                                                      theta,
                                                      pred_data,
                                                      samp_size,
                                                      num_responses[item],
                                                      num_quad)
      }
    }

    # Obtain E-table.
    for(i in 1:samp_size) {

      # Get weights for computing adaptive or fixed-point quadrature.
      posterior <- dnorm(theta,
                         mean = alpha[i],
                         sd = sqrt(phi[i]))

      # For each individual i, loop over items (j) and compute posterior
      # probability of response pattern.
      for(j in 1:num_items) {

        if(is.na(item_data[i,j])) next
        x <- item_data[i,j]

        if(item_type[item] == "cfa") { # Continuous responses.
          posterior <- posterior*itemtrace[[j]][i,]
        } else if(num_responses[j] == 2) { # Binary responses.
          if(x == 1) {
            posterior <- posterior*(1-itemtrace[[j]][i,])
          } else {
            posterior <- posterior*itemtrace[[j]][i,]
          }
        } else { # Ordered categorical responses.
          if(x == 1) {
            posterior <- posterior*(1-itemtrace[[j]][[1]][i,])
          } else if(x == num_responses[j]) {
            posterior <- posterior*itemtrace[[j]][[num_responses[j]-1]][i,]
          } else {
            posterior <- posterior*(itemtrace[[j]][[x-1]][i,]-
                                        itemtrace[[j]][[x]][i,])
          }
        }
      }

      # Normalize posterior.
      marginal <- sum(posterior, na.rm = TRUE)
      observed_ll <- observed_ll + log(marginal)
      if(marginal == 0) marginal <- 1
      etable[i,] <- posterior/marginal

      # EAPs.
      if(get_eap) {
        eap_scores[which(!NA_cases)[i]] <- sum(posterior*theta)/marginal
        eap_sd[which(!NA_cases)[i]] <-
          sqrt(sum(posterior*(theta-eap_scores[which(!NA_cases)[i]])**2)/marginal)
      }

    }

    # E-table matrix to be used in Q function and (possibly adaptive)
    # theta values.
    return(list(etable=etable,eap_scores=eap_scores,eap_sd=eap_sd,
                theta=theta,observed_ll=observed_ll))


  }
