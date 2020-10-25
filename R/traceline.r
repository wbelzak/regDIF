#' Binary item tracelines.
#'
#' @param p_active Vector of item parameters.
#' @param theta Matrix of adaptive theta values.
#' @param predictor.data Matrix or dataframe of DIF and/or impact predictors.
#' @param samp_size Sample size in dataset.
#' @param num_quadpts Number of quadrature points used for approximating the
#' latent variable.
#'
#' @keywords internal
#'
bernoulli_traceline_pts <-
  function(p_active,
           theta,
           predictor.data,
           samp_size,
           num_quadpts) {

    c0_parms <- grepl("c0",names(p_active),fixed=T)
    c1_parms <- grepl("c1",names(p_active),fixed=T)
    a0_parms <- grepl("a0",names(p_active),fixed=T)
    a1_parms <- grepl("a1",names(p_active),fixed=T)

    traceline0 <-
      apply(theta,
            2,
            function(x) {
              1 - 1 / (1 + exp(-((p_active[c0_parms] +
                                  predictor.data %*% p_active[c1_parms]) +
                                   (p_active[a0_parms] +
                                      predictor.data %*% p_active[a1_parms])*x)
                               )
                       )
              })
    traceline1 <- -traceline0 + 1

    return(list(traceline0,traceline1))

  }

#' Ordinal item tracelines (graded response).
#'
#' @param p_active Vector of item parameters.
#' @param theta Matrix of adaptive theta values.
#' @param predictor.data Matrix or dataframe of DIF and/or impact predictors.
#' @param samp_size Sample size in dataset.
#' @param num_responses_item Number of responses for item.
#' @param num_quadpts Number of quadrature points used for approximating the
#' latent variable.
#'
#' @keywords internal
#'
categorical_traceline_pts <-
  function(p_active,
           theta,
           predictor.data,
           samp_size,
           num_responses_item,
           num_quadpts) {

  # Space for category traceline (y = c category).
  traceline <- replicate(n=num_responses_item,
                         matrix(0,nrow=samp_size,ncol=num_quadpts),
                         simplify = F)

  c0_parms <- grepl("c0",names(p_active),fixed=T)
  c1_parms <- grepl("c1",names(p_active),fixed=T)
  a0_parms <- grepl("a0",names(p_active),fixed=T)
  a1_parms <- grepl("a1",names(p_active),fixed=T)

  # For item response 1.
  traceline[[1]] <-
    apply(theta,
          2,
          function(x) {
            1 - 1 / (1 + exp(-((p_active[c0_parms][1] +
                                  predictor.data %*% p_active[c1_parms]) +
                                 (p_active[a0_parms] +
                                    predictor.data %*% p_active[a1_parms])*x)
                             )
                     )
            })

  #for item response 2
  if(num_responses_item > 2) {
  traceline[[2]] <-
    apply(theta,
          2,
          function(x) {
            1 / (1 + exp(-((p_active[c0_parms][1] +
                              predictor.data %*% p_active[c1_parms]) +
                             (p_active[a0_parms] +
                                predictor.data %*% p_active[a1_parms])*x)
                         )
                 )
            }) -
    apply(theta,
          2,
          function(x) {
            1 / (1 + exp(-((p_active[c0_parms][1] -
                              p_active[c0_parms][2] +
                              predictor.data %*% p_active[c1_parms]) +
                             (p_active[a0_parms] +
                                predictor.data %*% p_active[a1_parms])*x)
                         )
                 )
            })

    # For item responses 3 to J-1.
    if(num_responses_item > 3) {
      for(thr in 3:(num_responses_item-1)) {
        traceline[[thr]] <-
          apply(theta,
                2,
                function(x) {
                  1 / (1 + exp(-((p_active[c0_parms][1] -
                                    p_active[c0_parms][thr-1] +
                                    predictor.data %*% p_active[c1_parms]) +
                                   (p_active[a0_parms] +
                                      predictor.data %*% p_active[a1_parms])*x)
                               )
                       )
                  }) -
          apply(theta,
                2,
                function(x) {
                  1 / (1 + exp(-((p_active[c0_parms][1] -
                                    p_active[c0_parms][thr] +
                                    predictor.data %*% p_active[c1_parms]) +
                                   (p_active[a0_parms] +
                                      predictor.data %*% p_active[a1_parms])*x)
                               )
                       )
                  })

      }

    }

  }

  # For item response J.
  traceline[[num_responses_item]] <-
    apply(theta,
          2,
          function(x) {
            1 / (1 + exp(-((p_active[c0_parms][1] -
                              p_active[c0_parms][num_responses_item-1] +
                              predictor.data %*% p_active[c1_parms]) +
                             (p_active[a0_parms] +
                                predictor.data %*% p_active[a1_parms])*x)
                         )
                 )
            })

  return(traceline)

  }

#' Ordinal tracelines (for derivatives).
#'
#' @param p_active Vector of item parameters.
#' @param theta Matrix of adaptive theta values.
#' @param predictor.data Matrix or dataframe of DIF and/or impact predictors.
#' @param samp_size Sample size in dataset.
#' @param num_responses_item Number of responses for item.
#' @param num_quadpts Number of quadrature points used for approximating the
#' latent variable.
#'
#' @keywords internal
#'
cumulative_traceline_pts <-
  function(p_active,
           theta,
           predictor.data,
           samp_size,
           num_responses_item,
           num_quadpts) {

  # Space for cumulative traceline (y >= c category).
  traceline <-
    replicate(n=(num_responses_item-1),
              matrix(0,nrow=samp_size,ncol=num_quadpts),
              simplify = F)

  c0_parms <- grepl("c0",names(p_active),fixed=T)
  c1_parms <- grepl("c1",names(p_active),fixed=T)
  a0_parms <- grepl("a0",names(p_active),fixed=T)
  a1_parms <- grepl("a1",names(p_active),fixed=T)

  # For item response 1.
  traceline[[1]] <-
    apply(theta,
          2,
          function(x) {
            1 / (1 + exp(-((p_active[c0_parms][1] +
                              predictor.data %*% p_active[c1_parms]) +
                             (p_active[a0_parms] +
                                predictor.data %*% p_active[a1_parms])*x)
                         )
                 )
            })

  # For item response 2 to J.
  if(num_responses_item > 2) {
    for(thr in 2:(num_responses_item-1)) {
      traceline[[thr]] <-
        apply(theta,
              2,
              function(x) {
                1 / (1 + exp(-((p_active[c0_parms][1] -
                                  p_active[c0_parms][thr] +
                                  predictor.data %*% p_active[c1_parms]) +
                                 (p_active[a0_parms] +
                                    predictor.data %*% p_active[a1_parms])*x)
                             )
                     )
                })
    }
  }

  return(traceline)

}

#' Continuous tracelines.
#'
#' @param p_active Vector of item parameters.
#' @param theta Matrix of adaptive theta values.
#' @param responses_item Vector of item responses.
#' @param predictor.data Matrix or dataframe of DIF and/or impact predictors.
#' @param samp_size Sample size in dataset.
#' @param num_quadpts Number of quadrature points used for approximating the
#' latent variable.
#'
#' @keywords internal
#'
gaussian_traceline_pts <-
  function(p_active,
           theta,
           responses_item,
           predictor.data,
           samp_size,
           num_quadpts) {

  c0_parms <- grepl("c0",names(p_active),fixed=T)
  c1_parms <- grepl("c1",names(p_active),fixed=T)
  a0_parms <- grepl("a0",names(p_active),fixed=T)
  a1_parms <- grepl("a1",names(p_active),fixed=T)
  s0_parms <- grepl("s0",names(p_active),fixed=T)
  s1_parms <- grepl("s1",names(p_active),fixed=T)

  mu <-
    apply(theta,
          2,
          function(x) {
            (p_active[c0_parms] +
               predictor.data %*% p_active[c1_parms]) +
              (p_active[a0_parms] +
                 predictor.data %*% p_active[a1_parms])*x
            })
  sigma <-
    sqrt(p_active[s0_parms][1]*exp(predictor.data %*% p_active[s1_parms]))

  traceline <- t(sapply(1:samp_size,
                        function(x) {
                          dnorm(responses_item[x],mu[x,],sigma[x])
                          }
                        ))

  return(list(traceline))

}



