#' Binary item tracelines.
#'
#' @param p_item Vector of item parameters.
#' @param theta Vector of theta values.
#' @param pred_data Matrix or dataframe of DIF and/or impact predictors.
#' @param samp_size Sample size in dataset.
#'
#' @return a \code{"matrix"} of probability values for Bernoulli item likelihood
#'
#' @keywords internal
#'
bernoulli_traceline_pts <-
  function(p_item,
           theta,
           pred_data,
           samp_size) {

    traceline <-
      vapply(theta,
             function(x) {
               1 / (1 + exp(
                 -((p_item[1] + pred_data %*% p_item[3:(2+ncol(pred_data))]) +
                   (p_item[2] + pred_data %*% p_item[(3+ncol(pred_data)):length(p_item)])*x)
                 ))
               }, numeric(samp_size))

    return(traceline)

  }

#' Binary item tracelines.
#'
#' @param p_item Vector of item parameters.
#' @param theta Vector of theta values.
#' @param pred_data Matrix or dataframe of DIF and/or impact predictors.
#' @param samp_size Sample size in dataset.
#'
#' @return a \code{"matrix"} of probability values for Bernoulli item likelihood
#'
#' @keywords internal
#'
bernoulli_traceline_pts2 <-
  function(p_item,
           theta,
           pred_data,
           samp_size) {

    traceline <-
      vapply(theta,
             function(x) {
               1 / (1 + exp(-(
                 (p_item[1] + pred_data %*% p_item[3:(2+ncol(pred_data))]) +
                 (p_item[2] + pred_data %*% p_item[(3+ncol(pred_data)):length(p_item)])*x
                 )))
             }, numeric(samp_size))

    return(colMeans(traceline))

  }

#' Binary item tracelines for proxy scores.
#'
#' @param p_item Vector of item parameters.
#' @param prox_data Vector of observed proxy scores.
#' @param pred_data Matrix or dataframe of DIF and/or impact predictors.
#'
#' @return a \code{"matrix"} of probability values for Bernoulli item likelihood using observed
#' proxy scores
#'
#' @keywords internal
#'
bernoulli_traceline_pts_proxy <-
  function(p_item,
           prox_data,
           pred_data) {

    traceline <-
      1 / (1 + exp(
        -((p_item[1] + pred_data %*% p_item[3:(2+ncol(pred_data))]) +
            (p_item[2] + pred_data %*% p_item[(3+ncol(pred_data)):length(p_item)])*prox_data)))

    return(traceline)

  }

#' Ordinal tracelines.
#'
#' @param p_item Vector of item parameters.
#' @param theta Vector of theta values.
#' @param pred_data Matrix or dataframe of DIF and/or impact predictors.
#' @param samp_size Sample size in dataset.
#' @param num_responses_item Number of responses for item.
#' @param num_quad Number of quadrature points used for approximating the
#' latent variable.
#'
#' @return a \code{"matrix"} of probability values for categorical (cumulative) item likelihood
#'
#' @keywords internal
#'
cumulative_traceline_pts <-
  function(p_item,
           theta,
           pred_data,
           samp_size,
           num_responses_item,
           num_quad) {

  # Space for cumulative traceline (y >= c category).
  traceline <-
    lapply(1:(num_responses_item-1), function(x) {
              matrix(0,nrow=samp_size,ncol=num_quad)})

  c0_parms <- grepl("c0",names(p_item),fixed=T)
  c1_parms <- grepl("c1",names(p_item),fixed=T)
  a0_parms <- grepl("a0",names(p_item),fixed=T)
  a1_parms <- grepl("a1",names(p_item),fixed=T)

  # For item response 1.
  traceline[[1]] <-
    vapply(theta,
          function(x) {
            1 / (1 + exp(-((p_item[c0_parms][1] +
                              pred_data %*% p_item[c1_parms]) +
                             (p_item[a0_parms] +
                                pred_data %*% p_item[a1_parms])*x)
                         )
                 )
            }, numeric(samp_size))

  # For item response 2 to J.
    for(thr in 2:(num_responses_item-1)) {
      traceline[[thr]] <-
        vapply(theta,
              function(x) {
                1 / (1 + exp(-((p_item[c0_parms][1] -
                                  p_item[c0_parms][thr] +
                                  pred_data %*% p_item[c1_parms]) +
                                 (p_item[a0_parms] +
                                    pred_data %*% p_item[a1_parms])*x)
                             )
                     )
                }, numeric(samp_size))
    }

  return(traceline)

}

#' Ordinal tracelines using proxy data.
#'
#' @param p_item Vector of item parameters.
#' @param prox_data Vector of observed proxy scores.
#' @param pred_data Matrix or dataframe of DIF and/or impact predictors.
#' @param samp_size Sample size in dataset.
#' @param num_responses_item Number of responses for item.
#'
#' @return a \code{"matrix"} of probability values for categorical (cumulative) item likelihood
#'
#' @keywords internal
#'
cumulative_traceline_pts_proxy <-
  function(p_item,
           prox_data,
           pred_data,
           samp_size,
           num_responses_item) {

    # Space for cumulative traceline (y >= c category).
    traceline <- lapply(1:(num_responses_item-1), function(x) {
                  matrix(0,nrow=samp_size,ncol=1)})

    c0_parms <- grepl("c0",names(p_item),fixed=T)
    c1_parms <- grepl("c1",names(p_item),fixed=T)
    a0_parms <- grepl("a0",names(p_item),fixed=T)
    a1_parms <- grepl("a1",names(p_item),fixed=T)

    # For item response 1.
    traceline[[1]] <- 1 / (1 + exp(-((p_item[c0_parms][1] +
                                 pred_data %*% p_item[c1_parms]) +
                                (p_item[a0_parms] +
                                   pred_data %*% p_item[a1_parms])*prox_data))
                           )

    # For item response 2 to J.
    for(thr in 2:(num_responses_item-1)) {
      traceline[[thr]] <- 1 / (1 + exp(-((p_item[c0_parms][1] -
                                   p_item[c0_parms][thr] +
                                   pred_data %*% p_item[c1_parms]) +
                                  (p_item[a0_parms] +
                                     pred_data %*% p_item[a1_parms])*prox_data))
                               )
    }

    return(traceline)

  }


#' Continuous tracelines.
#'
#' @param p_item Vector of item parameters.
#' @param theta Vector of theta values.
#' @param responses_item Vector of item responses.
#' @param pred_data Matrix or data frame of DIF and/or impact predictors.
#' @param samp_size Sample size in data set.
#'
#' @return a \code{"matrix"} of probability values for Gaussian item likelihood
#'
#' @keywords internal
#'
gaussian_traceline_pts <-
  function(p_item,
           theta,
           responses_item,
           pred_data,
           samp_size) {

  c0_parms <- grepl("c0",names(p_item),fixed=T)
  c1_parms <- grepl("c1",names(p_item),fixed=T)
  a0_parms <- grepl("a0",names(p_item),fixed=T)
  a1_parms <- grepl("a1",names(p_item),fixed=T)
  s0_parms <- grepl("s0",names(p_item),fixed=T)
  s1_parms <- grepl("s1",names(p_item),fixed=T)

  mu <-
    vapply(theta,
          function(x) {
            (p_item[c0_parms] +
               pred_data %*% p_item[c1_parms]) +
              (p_item[a0_parms] +
                 pred_data %*% p_item[a1_parms])*x
            }, numeric(samp_size))
  sigma <-
    sqrt(p_item[s0_parms][1]*exp(pred_data %*% p_item[s1_parms]))

  traceline <- t(sapply(1:samp_size,
                        function(x) {
                          dnorm(responses_item[x],mu[x,],sigma[x])
                          }
                        ))

  return(traceline)

}

#' Continuous tracelines using proxy data.
#'
#' @param p_item Vector of item parameters.
#' @param prox_data Vector of observed proxy scores.
#' @param responses_item Vector of item responses.
#' @param pred_data Matrix or data frame of DIF and/or impact predictors.
#' @param samp_size Sample size in data set.
#'
#' @return a \code{"matrix"} of probability values for Gaussian item likelihood
#'
#' @keywords internal
#'
gaussian_traceline_pts_proxy <-
  function(p_item,
           prox_data,
           responses_item,
           pred_data,
           samp_size) {

    c0_parms <- grepl("c0",names(p_item),fixed=T)
    c1_parms <- grepl("c1",names(p_item),fixed=T)
    a0_parms <- grepl("a0",names(p_item),fixed=T)
    a1_parms <- grepl("a1",names(p_item),fixed=T)
    s0_parms <- grepl("s0",names(p_item),fixed=T)
    s1_parms <- grepl("s1",names(p_item),fixed=T)

    mu <- (p_item[c0_parms] +
                  pred_data %*% p_item[c1_parms]) +
                 (p_item[a0_parms] +
                    pred_data %*% p_item[a1_parms])*prox_data
    sigma <-
      sqrt(p_item[s0_parms][1]*exp(pred_data %*% p_item[s1_parms]))

    traceline <- t(sapply(1:samp_size,
                          function(x) {
                            dnorm(responses_item[x],mu[x,],sigma[x])
                          }
    ))

    return(traceline)

  }

