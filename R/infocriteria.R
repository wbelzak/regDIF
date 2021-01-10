#' Maximization step.
#'
#' @param etable Matrix of E-table values for item and impact equations.
#' @param p List of parameters.
#' @param item.data Matrix or dataframe of item responses.
#' @param pred.data Matrix or dataframe of DIF and/or impact predictors.
#' @param mean_predictors Possibly different matrix of predictors for the mean
#' impact equation.
#' @param var_predictors Possibly different matrix of predictors for the
#' variance impact equation.
#' @param gamma Numeric value indicating the gamma parameter in the MCP
#' function.
#' @param samp_size Sample size in data set.
#' @param num_responses Number of responses for each item.
#' @param num_items Number of items in data set.
#' @param num.quad Number of quadrature points used for approximating the
#' latent variable.
#'
#' @keywords internal
#'
information_criteria <-
  function(etable,
           p,
           item.data,
           pred.data,
           mean_predictors,
           var_predictors,
           gamma,
           samp_size,
           num_responses,
           num_items,
           num.quad) {

  # Update theta and etable.
  theta <- etable$theta
  etable <- etable$etable

  ll_dif <- 0
  for (item in 1:num_items) {

    # Obtain E-tables for each response category.
    etable_item <- replicate(n=num_responses[item],
                        etable,
                        simplify = F)
    for(resp in 1:num_responses[item]) {
      etable_item[[resp]][which(
        !(item.data[,item] == resp)), ] <- 0
    }

    #compute negative log-likelihood values
    if(num_responses[item] == 1) {
      itemtrace <- gaussian_traceline_pts(p[[item]],
                                          theta,
                                          item.data[,item],
                                          pred.data,
                                          samp_size)
      ll_dif_item <- -1*sum(etable[[1]]*log(itemtrace[[1]]), na.rm = TRUE)

    } else if (num_responses[item] == 2) {
      itemtrace <- bernoulli_traceline_cpp(p[[item]],
                                           theta,
                                           pred.data,
                                           samp_size,
                                           num.quad)
      ll_dif_item <- -1*(sum(etable_item[[1]]*log(1-itemtrace), na.rm = TRUE) +
                           sum(etable_item[[2]]*log(itemtrace), na.rm = TRUE))

    } else if (num_responses[item] > 2){
      itemtrace <- cumulative_traceline_pts(p[[item]],
                                             theta,
                                             pred.data,
                                             samp_size,
                                             num_responses[item],
                                             num.quad)
      ll_dif_item <- 0
      for(resp in 1:num_responses[item]){
      if(all(itemtrace[[resp]] == 0)){
        log_itemtrace_cat <- 0
      } else{
        if(resp == 1) {
          log_itemtrace_cat <- log(1-itemtrace[[resp]])
        } else if(resp == num_responses[item]) {
          log_itemtrace_cat <- log(itemtrace[[resp]])
        } else {
          log_itemtrace_cat <- log(itemtrace[[resp-1]]-itemtrace[[resp]])
        }
      }
        log_itemtrace_cat[is.infinite(log_itemtrace_cat)] <- NA
        ll_dif_item <- ll_dif_item +
          -1*sum(etable[[resp]]*log_itemtrace_cat, na.rm = TRUE)
      }
    }

    ll_dif <- ll_dif + ll_dif_item
  }

  # Obtain likelihood value for latent variable model
  alpha <- mean_predictors %*% p[[num_items+1]]
  phi <- exp(var_predictors %*% p[[num_items+2]])
  prior_scores <- t(sapply(1:samp_size,
                           function(x) {
                             dnorm(theta,
                                   mean = alpha[x],
                                   sd = sqrt(phi[x]))
                             }))
  ll_impact <- -1*sum(etable*log(prior_scores))

  # Remove DIF parameters that equal zero from information crit calculation.
  p2 <- unlist(p)
  p2_base <- p2[c(grep("c0",names(p2)),grep("a0",names(p2)))]
  p2_cov <- p2[c(grep("c1",names(p2)),grep("a1",names(p2)))]
  p2_cov <- p2_cov[p2_cov != 0]
  p2 <- c(p2_base,p2_cov,p[[num_items+1]],p[[num_items+2]])
  ll <- ll_impact + ll_dif

  # Compute AIC and BIC.
  aic <- 2*length(p2) + 2*ll
  bic <- log(samp_size)*length(p2) + 2*ll

  return(c(aic,bic))

}
