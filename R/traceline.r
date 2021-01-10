#' Binary item tracelines.
#'
#' @param p_active Vector of item parameters.
#' @param theta Vector of theta values.
#' @param pred.data Matrix or dataframe of DIF and/or impact predictors.
#' @param samp_size Sample size in dataset.
#'
#' @keywords internal
#'
bernoulli_traceline_pts <-
  function(p_active,
           theta,
           pred.data,
           samp_size) {

    c0_parms <- grepl("c0",names(p_active),fixed=T)
    c1_parms <- grepl("c1",names(p_active),fixed=T)
    a0_parms <- grepl("a0",names(p_active),fixed=T)
    a1_parms <- grepl("a1",names(p_active),fixed=T)

    traceline <- vapply(theta,
                        function(x) {1 / (1 + exp(-((p_active[c0_parms] +
                                  pred.data %*% p_active[c1_parms]) +
                                   (p_active[a0_parms] +
                                      pred.data %*% p_active[a1_parms])*x)
                               )
                       )
              },numeric(samp_size))

    return(traceline)

  }

#' Ordinal tracelines (for derivatives).
#'
#' @param p_active Vector of item parameters.
#' @param theta Vector of theta values.
#' @param pred.data Matrix or dataframe of DIF and/or impact predictors.
#' @param samp_size Sample size in dataset.
#' @param num_responses_item Number of responses for item.
#' @param num.quad Number of quadrature points used for approximating the
#' latent variable.
#'
#' @keywords internal
#'
cumulative_traceline_pts <-
  function(p_active,
           theta,
           pred.data,
           samp_size,
           num_responses_item,
           num.quad) {

  # Space for cumulative traceline (y >= c category).
  traceline <-
    replicate(n=(num_responses_item-1),
              matrix(0,nrow=samp_size,ncol=num.quad),
              simplify = F)

  c0_parms <- grepl("c0",names(p_active),fixed=T)
  c1_parms <- grepl("c1",names(p_active),fixed=T)
  a0_parms <- grepl("a0",names(p_active),fixed=T)
  a1_parms <- grepl("a1",names(p_active),fixed=T)

  # For item response 1.
  traceline[[1]] <-
    vapply(theta,
          function(x) {
            1 / (1 + exp(-((p_active[c0_parms][1] +
                              pred.data %*% p_active[c1_parms]) +
                             (p_active[a0_parms] +
                                pred.data %*% p_active[a1_parms])*x)
                         )
                 )
            }, numeric(samp_size))

  # For item response 2 to J.
    for(thr in 2:(num_responses_item-1)) {
      traceline[[thr]] <-
        vapply(theta,
              function(x) {
                1 / (1 + exp(-((p_active[c0_parms][1] -
                                  p_active[c0_parms][thr] +
                                  pred.data %*% p_active[c1_parms]) +
                                 (p_active[a0_parms] +
                                    pred.data %*% p_active[a1_parms])*x)
                             )
                     )
                }, numeric(samp_size))
    }

  return(traceline)

}

#' Continuous tracelines.
#'
#' @param p_active Vector of item parameters.
#' @param theta Vector of theta values.
#' @param responses_item Vector of item responses.
#' @param pred.data Matrix or dataframe of DIF and/or impact predictors.
#' @param samp_size Sample size in dataset.
#'
#' @keywords internal
#'
gaussian_traceline_pts <-
  function(p_active,
           theta,
           responses_item,
           pred.data,
           samp_size) {

  c0_parms <- grepl("c0",names(p_active),fixed=T)
  c1_parms <- grepl("c1",names(p_active),fixed=T)
  a0_parms <- grepl("a0",names(p_active),fixed=T)
  a1_parms <- grepl("a1",names(p_active),fixed=T)
  s0_parms <- grepl("s0",names(p_active),fixed=T)
  s1_parms <- grepl("s1",names(p_active),fixed=T)

  mu <-
    vapply(theta,
          function(x) {
            (p_active[c0_parms] +
               pred.data %*% p_active[c1_parms]) +
              (p_active[a0_parms] +
                 pred.data %*% p_active[a1_parms])*x
            }, numeric(samp_size))
  sigma <-
    sqrt(p_active[s0_parms][1]*exp(pred.data %*% p_active[s1_parms]))

  traceline <- t(sapply(1:samp_size,
                        function(x) {
                          dnorm(responses_item[x],mu[x,],sigma[x])
                          }
                        ))

  return(list(traceline))

}



