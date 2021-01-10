#' Partial derivatives for mean impact equation.
#'
#' @param p_active Vector of item parameters.
#' @param etable_all E-table for impact.
#' @param theta Matrix of adaptive theta values.
#' @param mean_predictors Possibly different matrix of predictors for the mean
#' impact equation.
#' @param var_predictors Possibly different matrix of predictors for the
#' variance impact equation.
#' @param cov Covariate being maximized.
#' @param samp_size Sample size in data set.
#' @param num_items Number of items in data set.
#' @param num.quad Number of quadrature points used for approximating the
#' latent variable.
#'
#' @keywords internal
#'
d_alpha <-
  function(p_impact,
           etable,
           theta,
           mean_predictors,
           var_predictors,
           cov,
           samp_size,
           num_items,
           num.quad) {

  # Get latent mean and variance vectors.
  alpha <- mean_predictors %*% p_impact[grep("g",names(p_impact),fixed=T)]
  phi <- exp(var_predictors %*% p_impact[grep("b",names(p_impact),fixed=T)])

  eta_d <- matrix(rep(mean_predictors[,cov+1], num.quad),
                  ncol = num.quad,
                  nrow = samp_size)

  d1_trace <- vapply(1:num.quad,
                       function(x) {
                         eta_d[,x]/phi*(theta[x]-alpha)
                         },numeric(samp_size))
  d2_trace <- vapply(1:num.quad,
                       function(x) {
                         -eta_d[,x]**2/phi
                         },numeric(samp_size))

  d1 <- sum(etable*d1_trace, na.rm = TRUE)
  d2 <- sum(etable*d2_trace, na.rm = TRUE)

  dlist <- list(d1,d2)

}

#' Partial derivatives for mean impact equation.
#'
#' @param p_active Vector of item parameters.
#' @param etable_all E-table for impact.
#' @param theta Matrix of adaptive theta values.
#' @param mean_predictors Possibly different matrix of predictors for the mean
#' impact equation.
#' @param var_predictors Possibly different matrix of predictors for the
#' variance impact equation.
#' @param cov Covariate being maximized.
#' @param samp_size Sample size in dataset.
#' @param num_items Number of items in dataset.
#' @param num.quad Number of quadrature points used for approximating the
#' latent variable.
#'
#' @keywords internal
#'
d_phi <-
  function(p_impact,
           etable,
           theta,
           mean_predictors,
           var_predictors,
           cov,
           samp_size,
           num_items,
           num.quad) {

  # Get latent mean and variance vectors
  alpha <- mean_predictors %*% p_impact[grep("g",names(p_impact),fixed=T)]
  phi <- exp(var_predictors %*% p_impact[grep("b",names(p_impact),fixed=T)])

  eta_d1 <- .5*sqrt(phi)*var_predictors[,cov+1]
  eta_d2 <- .25*sqrt(phi)*var_predictors[,cov+1]**2

  d1_trace <- vapply(1:num.quad,
                       function(x) {
                         eta_d1*((theta[x]-alpha)**2/phi**(3/2) -
                                      1/sqrt(phi))
                         },numeric(samp_size))
  d2_trace <- vapply(1:num.quad,
                       function(x) {
                         -2*eta_d2*(phi**(-3/2)*(theta[x]-alpha)**2)
                         },numeric(samp_size))

  d1 <- sum(etable*d1_trace, na.rm = TRUE)
  d2 <- sum(etable*d2_trace, na.rm = TRUE)

  dlist <- list(d1,d2)

}

#' Partial derivatives for binary items.
#'
#' @param parm Item parameter being maximized.
#' @param p_item Vector of item parameters.
#' @param etable_item E-table for item.
#' @param theta Matrix of adaptive theta values.
#' @param pred.data Matrix or dataframe of DIF and/or impact predictors.
#' @param cov Covariate being maximized.
#' @param samp_size Sample size in dataset.
#' @param num_items Number of items in dataset.
#' @param num.quad Number of quadrature points used for approximating the
#' latent variable.
#'
#' @keywords internal
#'
d_bernoulli <-
  function(parm,
           p_item,
           etable_item,
           theta,
           pred.data,
           cov,
           samp_size,
           num_items,
           num.quad) {

  if(parm == "c0"){
    eta_d <- matrix(1, nrow = samp_size, ncol = num.quad)
  } else if(parm == "a0"){
    eta_d <- t(replicate(n=samp_size, theta))
  } else if(parm == "c1"){
    eta_d <- matrix(rep(pred.data[,cov+1], num.quad),
                    ncol = num.quad,
                    nrow = samp_size)
  } else if(parm == "a1"){
    eta_d <- matrix(rep(pred.data[,cov+1], num.quad),
                    ncol = num.quad,
                    nrow = samp_size)*t(replicate(n=samp_size, theta))
  }

  traceline <- bernoulli_traceline_pts(p_item,
                                       theta,
                                       pred.data,
                                       samp_size)

  d1 <- sum(eta_d*traceline*(etable_item[[2]]/traceline -
                               etable_item[[2]] -
                               etable_item[[1]]),
            na.rm = TRUE)
  d2 <- sum(eta_d**2*(-traceline + traceline**2)*(etable_item[[1]] +
                                                    etable_item[[2]]),
            na.rm = TRUE)

  dlist <- list(d1,d2)

  }

#' Partial derivatives for ordinal items.
#'
#' @param parm Item parameter being maximized.
#' @param p_item Vector of item parameters.
#' @param etable_item E-table for impact.
#' @param theta Matrix of adaptive theta values.
#' @param pred.data Matrix or dataframe of DIF and/or impact predictors.
#' @param thr Threshold value being maximized.
#' @param cov Covariate being maximized.
#' @param samp_size Sample size in dataset.
#' @param num_items Number of items in dataset.
#' @param num.quad Number of quadrature points used for approximating the
#' latent variable.
#'
#' @keywords internal
#'
d_categorical <-
  function(parm,
           p_item,
           etable_item,
           theta,
           pred.data,
           thr,
           cov,
           samp_size,
           num_responses_item,
           num_items,
           num.quad) {

    if(parm == "c0"){
      eta_d <- matrix(1, nrow = samp_size, ncol = num.quad)
    } else if(parm == "a0"){
      eta_d <- t(replicate(n=samp_size, theta))
    } else if(parm == "c1"){
      eta_d <- matrix(rep(pred.data[,cov], num.quad),
                      ncol = num.quad,
                      nrow = samp_size)
    } else if(parm == "a1"){
      eta_d <- matrix(rep(pred.data[,cov], num.quad),
                      ncol = num.quad,
                      nrow = samp_size)*t(replicate(n=samp_size, theta))
    }

    cum_traceline <- cumulative_traceline_pts(p_item,
                                              theta,
                                              pred.data,
                                              samp_size,
                                              num_responses_item,
                                              num.quad)

    # Non-threshold derivatives.
    if(thr < 0){
      d1 <- eta_d*(-etable_item[[1]]*cum_traceline[[1]] +
                     etable_item[[num_responses_item]]*(
                       1 - cum_traceline[[num_responses_item-1]]))
      d2 <- eta_d**2*(-etable_item[[1]]*(cum_traceline[[1]]*(
        1-cum_traceline[[1]])) +
          etable_item[[num_responses_item]]*(
            -cum_traceline[[num_responses_item-1]]*(
              1-cum_traceline[[num_responses_item-1]])))

      for(i in 2:(num_responses_item-1)){

        # Skip intermediate derivative calculations for constrained theshold.
        d1 <- d1 + eta_d*etable_item[[i]]*((1-cum_traceline[[i]]) -
                                             cum_traceline[[i-1]])
        d2 <- d2 + eta_d**2*etable_item[[i]]*(cum_traceline[[i-1]]**2 +
                                                cum_traceline[[i]]**2 -
                                                cum_traceline[[i-1]] -
                                                cum_traceline[[i]])
      }


      d1 <- sum(d1, na.rm = TRUE)
      d2 <- sum(d2, na.rm = TRUE)

      # Threshold derivatives.
    } else {
      if(thr < (num_responses_item-1)) {
        cat_traceline <- (cum_traceline[[thr]] - cum_traceline[[thr+1]])
      } else{
        cat_traceline <- cum_traceline[[thr]]
      }
      d1 <-
        sum(-etable_item[[thr]]*cum_traceline[[thr]]*(
          1 - cum_traceline[[thr]]
        ) / (cum_traceline[[thr-1]] - cum_traceline[[thr]]), na.rm = TRUE) +
        sum(etable_item[[thr+1]]*cum_traceline[[thr]]*(
          1 - cum_traceline[[thr]]
        ) / cat_traceline, na.rm = TRUE)
      d2 <- sum(etable_item[[thr]] / (cum_traceline[[thr-1]] -
                                        cum_traceline[[thr]])*(
                                          cum_traceline[[thr]]*(1 - cum_traceline[[thr]])**2 -
                                            cum_traceline[[thr]]**2*(1 - cum_traceline[[thr]]) +
                                            cum_traceline[[thr]]**2*(1 - cum_traceline[[thr]])**2 /
                                            (cum_traceline[[thr-1]] - cum_traceline[[thr]])
                                        ), na.rm = TRUE) -
        sum(etable_item[[thr+1]] / cat_traceline*(
                                        cum_traceline[[thr]]*(1 - cum_traceline[[thr]])**2 -
                                          cum_traceline[[thr]]**2*(1-cum_traceline[[thr]]) -
                                          cum_traceline[[thr]]**2*(1-cum_traceline[[thr]])**2 /
                                          cat_traceline
                                      ), na.rm = TRUE)

    }

    dlist <- list(d1,d2)

  }


#' Partial derivatives for mean parameter of continuous items.
#'
#' @param parm Item parameter being maximized.
#' @param p_item Vector of item parameters.
#' @param etable_item E-table for impact.
#' @param theta Matrix of adaptive theta values.
#' @param responses_item Vector of item responses.
#' @param pred.data Matrix or dataframe of DIF and/or impact predictors.
#' @param cov Covariate being maximized.
#' @param samp_size Sample size in dataset.
#' @param num_items Number of items in dataset.
#' @param num.quad Number of quadrature points used for approximating the
#' latent variable.
#'
#' @keywords internal
#'
d_mu_gaussian <-
  function(parm,
           p_item,
           etable_item,
           theta,
           responses_item,
           pred.data,
           cov,
           samp_size,
           num_items,
           num.quad) {

  if(parm == "c0"){
    eta_d <- matrix(1, nrow = samp_size, ncol = num.quad)
  } else if(parm == "a0"){
    eta_d <- t(replicate(n=samp_size, theta))
  } else if(parm == "c1"){
    eta_d <- matrix(rep(pred.data[,cov], num.quad),
                    ncol = num.quad,
                    nrow = samp_size)
  } else if(parm == "a1"){
    eta_d <- matrix(rep(pred.data[,cov], num.quad),
                    ncol = num.quad,
                    nrow = samp_size)*t(replicate(n=samp_size, theta))
  }


  # Get latent mean and variance vectors.
  mu <- vapply(theta,
              function(x) {
                (p_item[grep("c0",names(p_item),fixed=T)] +
                   pred.data %*%
                   p_item[grep("c1",names(p_item),fixed=T)]) +
                  (p_item[grep("a0",names(p_item),fixed=T)] +
                     pred.data %*%
                     p_item[grep("a1",names(p_item),fixed=T)])*x
                },numeric(samp_size))
  sigma <- sqrt(p_item[grep("s0",names(p_item))][1]*exp(
    pred.data %*% p_item[grep("s1",names(p_item))]
    ))


  d1_trace <- t(vapply(1:samp_size,
                       function(x) {
                         eta_d[x,]/sigma[x]**2*(responses_item[x] - mu[x,])
                         },numeric(samp_size)))
  d2_trace <- t(vapply(1:samp_size,
                       function(x) {
                         -eta_d[x,]**2 / sigma[x]**2
                         },numeric(samp_size)))

  d1 <- sum(etable_item[[1]]*d1_trace, na.rm = TRUE)
  d2 <- sum(etable_item[[1]]*d2_trace, na.rm = TRUE)

  dlist <- list(d1,d2)

}

#' Partial derivatives for variance parameter of continuous items.
#'
#' @param parm Item parameter being maximized.
#' @param p_item Vector of item parameters.
#' @param etable_item E-table for impact.
#' @param theta Matrix of adaptive theta values.
#' @param responses_item Vector of item responses.
#' @param pred.data Matrix or dataframe of DIF and/or impact predictors.
#' @param cov Covariate being maximized.
#' @param samp_size Sample size in dataset.
#' @param num_items Number of items in dataset.
#' @param num.quad Number of quadrature points used for approximating the
#' latent variable.
#'
#' @keywords internal
#'
d_sigma_gaussian <-
  function(parm,
           p_item,
           etable_item,
           theta,
           responses_item,
           pred.data,
           cov,
           samp_size,
           num_items,
           num.quad) {

  sigma <- sqrt(p_item[grep("s0",names(p_item))][1]*exp(
    pred.data %*% p_item[grep("s1",names(p_item))]))
  mu <- vapply(theta,
              function(x) {
                (p_item[grep("c0",names(p_item),fixed=T)] +
                   pred.data %*%
                   p_item[grep("c1",names(p_item),fixed=T)]) +
                  (p_item[grep("a0",names(p_item),fixed=T)] +
                     pred.data %*%
                     p_item[grep("a1",names(p_item),fixed=T)])*x
                },numeric(samp_size))

  if(parm == "s0") {
    eta_d1 <- vapply(1:samp_size,
                     function(x) {
                       exp(pred.data[x,] %*%
                             p_item[grep("s1",names(p_item))]) / (2*sigma[x])
                       },numeric(samp_size))
    eta_d2 <- vapply(1:samp_size,
                     function(x) {
                       -exp(pred.data[x,] %*%
                              p_item[grep("s1",names(p_item))])**2 /
                         (4*sigma[x]**3)
                       },numeric(samp_size))
  } else if(parm == "s1") {
    eta_d1 <- vapply(1:samp_size,
                     function(x) {
                       sigma[x]*pred.data[x,cov] / 2
                       },numeric(samp_size))
    eta_d2 <- vapply(1:samp_size,
                     function(x) {
                       sigma[x]*pred.data[x,cov]**2 / 4
                       },numeric(samp_size))
  }


  d1_trace <- t(vapply(1:samp_size,
                       function(x) {
                         eta_d1[x]*((responses_item[x]-mu[x,])**2 /
                                      sigma[x]**3 -
                                      1/sigma[x])
                         },numeric(samp_size)))

  if(parm == "s0") {
    d2_trace <- t(vapply(1:samp_size,
                         function(x) {
                           eta_d1[x]**2*(1 / sigma[x]**2 -
                                           3*(responses_item[x] - mu[x,])**2 /
                                           sigma[x]**4) +
                             eta_d2[x]*((responses_item[x] - mu[x,])**2 /
                                          sigma[x]**3 - 1/sigma[x])
                           },numeric(samp_size)))
  } else if(parm == "s1") {
    d2_trace <- t(vapply(1:samp_size,
                         function(x) {
                           -2*eta_d2[x]*(sigma[x]**(-3)*(responses_item[x] -
                                                           mu[x,])**2)
                           },numeric(samp_size)))
  }

  d1 <- sum(etable_item[[1]]*d1_trace, na.rm = TRUE)
  d2 <- sum(etable_item[[1]]*d2_trace, na.rm = TRUE)

  dlist <- list(d1,d2)

}



