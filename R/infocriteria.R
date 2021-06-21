#' Maximization step.
#'
#' @param eout E-table output.
#' @param p List of parameters.
#' @param item_data Matrix or data.frame of item responses.
#' @param pred_data Matrix or data.frame of DIF and/or impact predictors.
#' @param prox_data Vector of observed proxy scores.
#' @param mean_predictors Possibly different matrix of predictors for the mean
#' impact equation.
#' @param var_predictors Possibly different matrix of predictors for the
#' variance impact equation.
#' @param gamma Numeric value indicating the gamma parameter in the MCP
#' function.
#' @param samp_size Sample size in data set.
#' @param num_responses Number of responses for each item.
#' @param num_items Number of items in data set.
#' @param num_quad Number of quadrature points used for approximating the
#' latent variable.
#'
#' @return a \code{"list"} of information criteria to use for model selection
#'
#' @keywords internal
#'
information_criteria <-
  function(eout,
           p,
           item_data,
           pred_data,
           prox_data,
           mean_predictors,
           var_predictors,
           gamma,
           samp_size,
           num_responses,
           num_items,
           num_quad) {

  # Update theta and etable.
  if(!is.null(eout)) {
    theta <- eout$theta
    etable <- eout$etable
    theta_mat <- t(matrix(theta,
                          ncol=samp_size,
                          nrow=num_quad))
  }

  complete_ll_dif <- 0
  observed_ll_dif <- 0
  for (item in 1:num_items) {

    if(!is.null(eout)) {
      # Obtain E-tables for each response category.
      etable_item <- lapply(1:num_responses[item], function(x) etable)

      for(resp in 1:num_responses[item]) {
        etable_item[[resp]][which(
          !(item_data[,item] == resp)), ] <- 0
      }
    } else {
      # Obtain item data for each response category.
      item_data_resp <- lapply(1:num_responses[item], function(x) item_data)

      for(resp in 1:num_responses[item]) {

        item_data_resp[[resp]][!(item_data_resp[[resp]] == resp)] <- 0
        item_data_resp[[resp]][item_data_resp[[resp]] == resp] <- 1

      }


    }

    #compute negative log-likelihood values
    if(num_responses[item] == 1) {
      itemtrace <- gaussian_traceline_pts(p[[item]],
                                          theta,
                                          item_data[,item],
                                          pred_data,
                                          samp_size)
      complete_ll_dif_item <- sum(etable_item[[1]]*log(itemtrace[[1]]),
                                     na.rm = TRUE)
      observed_ll_dif_item <- sum(log(itemtrace[[1]]),
                                     na.rm = TRUE)

    } else if (num_responses[item] == 2) {

      if(!is.null(eout)) {
        itemtrace <- bernoulli_traceline_pts(p[[item]],
                                             theta,
                                             pred_data,
                                             samp_size)
        complete_ll_dif_item <- sum(etable_item[[1]]*log(1-itemtrace),
                                    na.rm = TRUE) + sum(etable_item[[2]]*log(itemtrace),
                                                        na.rm = TRUE)
        observed_ll_dif_item <- sum(theta_mat*log(1-itemtrace),
                                    na.rm = TRUE) + sum(theta_mat*log(itemtrace), na.rm = TRUE)

      } else {
        itemtrace <- bernoulli_traceline_pts_proxy(p[[item]],
                                                   prox_data,
                                                   pred_data)
        complete_ll_dif_item <- sum(item_data_resp[[1]][, item]*log(1-itemtrace), na.rm = TRUE) +
          sum(item_data_resp[[2]][, item]*log(itemtrace), na.rm = TRUE)
        observed_ll_dif_item <- sum(prox_data*log(1-itemtrace),
                                    na.rm = TRUE) + sum(prox_data*log(itemtrace), na.rm = TRUE)
      }


    } else if (num_responses[item] > 2){
      itemtrace <- cumulative_traceline_pts(p[[item]],
                                             theta,
                                             pred_data,
                                             samp_size,
                                             num_responses[item],
                                             num_quad)
      complete_ll_dif_item <- 0
      observed_ll_dif_item <- 0
      for(resp in 1:num_responses[item]){
        if(resp < num_responses[item] && all(itemtrace[[resp]] == 0)){
          log_itemtrace_cat <- 0
        } else{
          if(resp == 1) {
            log_itemtrace_cat <- log(1-itemtrace[[resp]])
          } else if(resp == num_responses[item]) {
            log_itemtrace_cat <- log(itemtrace[[resp-1]])
          } else {
            log_itemtrace_cat <- log(itemtrace[[resp-1]]-itemtrace[[resp]])
          }
        }
        log_itemtrace_cat[is.infinite(log_itemtrace_cat)] <- NA
        complete_ll_dif_item <- complete_ll_dif_item +
          sum(etable_item[[resp]]*log_itemtrace_cat, na.rm = TRUE)
        observed_ll_dif_item <- observed_ll_dif_item +
          sum(log_itemtrace_cat, na.rm = TRUE)
      }
    }

    complete_ll_dif <- complete_ll_dif + complete_ll_dif_item
    observed_ll_dif <- observed_ll_dif + observed_ll_dif_item
  }

  # Obtain likelihood value for latent variable model
  alpha <- mean_predictors %*% p[[num_items+1]]
  phi <- exp(var_predictors %*% p[[num_items+2]])

  if(!is.null(eout)) {
    prior_scores <- t(sapply(1:samp_size,
                             function(x) {
                               dnorm(theta,
                                     mean = alpha[x],
                                     sd = sqrt(phi[x]))
                             }))
    complete_ll_impact <- sum(etable*log(prior_scores), na.rm = TRUE)
  } else {
    prior_scores <- dnorm(prox_data,
                          mean = alpha,
                          sd = sqrt(phi))
    complete_ll_impact <- sum(log(prior_scores), na.rm = TRUE)
  }

  observed_ll_impact <- sum(log(prior_scores), na.rm = TRUE)

  # Remove DIF parameters that equal zero from information crit calculation.
  p2 <- unlist(p)
  p2_base <- p2[c(grep("c0",names(p2)),grep("a0",names(p2)))]
  p2_cov <- p2[c(grep("c1",names(p2)),grep("a1",names(p2)))]
  p2_cov <- p2_cov[p2_cov != 0]
  p2 <- c(p2_base,p2_cov,p[[num_items+1]],p[[num_items+2]])
  complete_ll <- complete_ll_impact + complete_ll_dif
  observed_ll <- observed_ll_impact + observed_ll_dif

  # Compute AIC and BIC.
  aic <- 2*length(p2) - 2*complete_ll
  bic <- log(samp_size)*length(p2) - 2*complete_ll

  return(list(aic=aic,bic=bic,complete_ll=complete_ll,observed_ll=observed_ll))

}
