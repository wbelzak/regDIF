#' Maximization step.
#'
#' @param elist List of E-tables for item and impact equations, in addition to
#' theta values.
#' @param p List of parameters.
#' @param item.data Matrix or dataframe of item responses.
#' @param pred.data Matrix or dataframe of DIF and/or impact predictors.
#' @param mean_predictors Possibly different matrix of predictors for the mean
#' impact equation.
#' @param var_predictors Possibly different matrix of predictors for the
#' variance impact equation.
#' @param theta Matrix of adaptive theta values.
#' @param gamma Numeric value indicating the gamma parameter in the MCP
#' function.
#' @param samp_size Sample size in dataset.
#' @param num_responses Number of responses for each item.
#' @param num_items Number of items in dataset.
#' @param quad.pts Number of quadrature points used for approximating the
#' latent variable.
#'
#' @keywords internal
#'
information_criteria <-
  function(elist,
           p,
           item.data,
           pred.data,
           mean_predictors,
           var_predictors,
           theta,
           gamma,
           samp_size,
           num_responses,
           num_items,
           quad.pts) {

  ll_dif <- 0
  for (item in 1:num_items) {

    # Obtain E-tables for each response category.
    etable <- replicate(n=num_responses[item],
                        elist[[1]][[item]],
                        simplify = F)
    for(resp in 1:num_responses[item]) {
      etable[[resp]][which(
        !(etable[[resp]][,ncol(etable[[resp]])] == resp)),] <- 0
    }
    etable <- lapply(etable, function(x) x[,1:quad.pts])

    #get item parameters
    p_item <- p[[item]]

    #compute negative log-likelihood values
    if(num_responses[item] == 1) {
      itemtrace <- gaussian_traceline_pts(p[[item]],
                                          theta,
                                          item.data[,item],
                                          pred.data,
                                          samp_size,
                                          quad.pts)
      ll_dif_item <- -1*sum(etable[[1]]*log(itemtrace[[1]]), na.rm = TRUE)

    } else if (num_responses[item] == 2) {
      itemtrace <- bernoulli_traceline_pts(p[[item]],
                                           theta,
                                           pred.data,
                                           samp_size,
                                           quad.pts)
      ll_dif_item <- -1*(sum(etable[[1]]*log(itemtrace[[1]]), na.rm = TRUE) +
                           sum(etable[[2]]*log(itemtrace[[2]]), na.rm = TRUE))

    } else if (num_responses[item] > 2){
      itemtrace <- categorical_traceline_pts(p[[item]],
                                             theta,
                                             pred.data,
                                             samp_size,
                                             num_responses[item],
                                             quad.pts)
      ll_dif_item <- 0
      for(resp in 1:num_responses[item]){
      if(all(itemtrace[[resp]] == 0)){
        log_itemtrace_cat <- 0
      } else{
        log_itemtrace_cat <- log(itemtrace[[resp]])
      }
        log_itemtrace_cat[is.infinite(log_itemtrace_cat)] <- NA
        ll_dif_item <- ll_dif_item +
          -1*sum(etable[[resp]]*log_itemtrace_cat, na.rm = TRUE)
      }
    }

    ll_dif <- ll_dif + ll_dif_item
  }

  # Obtain likelihood value for latent variable model
  p_impact <- c(p[[num_items+1]],p[[num_items+2]])
  alpha <- mean_predictors %*% p_impact[grep("g",names(p_impact),fixed=T)]
  phi <- exp(var_predictors %*% p_impact[grep("b",names(p_impact),fixed=T)])
  prior_scores <- t(sapply(1:samp_size,
                           function(x) {
                             dnorm(theta[x,],
                                   mean = alpha[x],
                                   sd = sqrt(phi[x]))
                             }))
  ll_impact <- -1*sum(elist[[2]]*log(prior_scores))

  # Remove DIF parameters that equal zero from information crit calculation.
  p2 <- unlist(p)
  p2_base <- p2[c(grep("c0",names(p2)),grep("a0",names(p2)))]
  p2_cov <- p2[c(grep("c1",names(p2)),grep("a1",names(p2)))]
  p2_cov <- p2_cov[p2_cov != 0]
  p2 <- c(p2_base,p2_cov,p_impact)
  ll <- ll_impact + ll_dif

  # Compute AIC and BIC.
  aic <- 2*length(p2) + 2*ll
  bic <- log(samp_size)*length(p2) + 2*ll

  return(c(aic,bic))

}
