#' Partial derivatives for mean impact equation.
#'
#' @param p_impact Vector of impact parameters.
#' @param etable E-table for impact.
#' @param theta Matrix of adaptive theta values.
#' @param mean_predictors Possibly different matrix of predictors for the mean
#' impact equation.
#' @param var_predictors Possibly different matrix of predictors for the
#' variance impact equation.
#' @param cov Covariate being maximized.
#' @param samp_size Sample size in data set.
#' @param num_items Number of items in data set.
#' @param num_quad Number of quadrature points used for approximating the
#' latent variable.
#'
#' @return a \code{"list"} of first and second partial derivatives for mean impact equation (to
#' use with coordinate descent and univariate Newton-Raphson)
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
           num_quad) {

  # Get latent mean and variance vectors.
  alpha <- mean_predictors %*% p_impact[grep("mean",names(p_impact),fixed=T)]
  phi <- exp(var_predictors %*% p_impact[grep("var",names(p_impact),fixed=T)])

  d1_trace <- vapply(1:num_quad,
                       function(x) {
                         mean_predictors[,cov]/phi*(theta[x]-alpha)
                         },numeric(samp_size))
  d2_trace <- vapply(1:num_quad,
                       function(x) {
                         -mean_predictors[,cov]**2/phi
                         },numeric(samp_size))

  d1 <- sum(etable*d1_trace, na.rm = TRUE)
  d2 <- sum(etable*d2_trace, na.rm = TRUE)

  dlist <- list(d1,d2)

}

#' Partial derivatives for mean impact equation.
#'
#' @param p_impact Vector of impact parameters.
#' @param etable E-table for impact.
#' @param theta Matrix of adaptive theta values.
#' @param mean_predictors Possibly different matrix of predictors for the mean
#' impact equation.
#' @param var_predictors Possibly different matrix of predictors for the
#' variance impact equation.
#' @param cov Covariate being maximized.
#' @param samp_size Sample size in dataset.
#' @param num_items Number of items in dataset.
#' @param num_quad Number of quadrature points used for approximating the
#' latent variable.
#'
#' @return a \code{"list"} of first and second partial derivatives for variance impact equation (to
#' use with coordinate descent and univariate Newton-Raphson)
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
           num_quad) {

  # Get latent mean and variance vectors
  alpha <- mean_predictors %*% p_impact[grep("mean",names(p_impact),fixed=T)]
  phi <- exp(var_predictors %*% p_impact[grep("var",names(p_impact),fixed=T)])

  eta_d1 <- .5*sqrt(phi)*var_predictors[,cov]
  eta_d2 <- .5*sqrt(phi)*var_predictors[,cov]**2

  d1_trace <- vapply(1:num_quad,
                       function(x) {
                         eta_d1*((theta[x]-alpha)**2/phi**(3/2) -
                                      1/sqrt(phi))
                         },numeric(samp_size))
  d2_trace <- vapply(1:num_quad,
                       function(x) {
                         -eta_d2*(phi**(-3/2)*(theta[x]-alpha)**2)
                         },numeric(samp_size))

  d1 <- sum(etable*d1_trace, na.rm = TRUE)
  d2 <- sum(etable*d2_trace, na.rm = TRUE)

  dlist <- list(d1,d2)

  }

#' Partial derivatives for mean and variance impact equation.
#'
#' @param p_mean Vector of mean impact parameters.
#' @param p_var Vector of variance impact parameters.
#' @param etable E-table for impact.
#' @param theta Matrix of adaptive theta values.
#' @param mean_predictors Possibly different matrix of predictors for the mean
#' impact equation.
#' @param var_predictors Possibly different matrix of predictors for the
#' variance impact equation.
#' @param samp_size Sample size in data set.
#' @param num_items Number of items in data set.
#' @param num_quad Number of quadrature points used for approximating the
#' latent variable.
#' @param num_predictors Number of predictors in dataset.
#'
#' @return a \code{"list"} of first and second partial derivatives for impact equation (to use
#' with multivariate Newton-Raphson)
#'
#' @keywords internal
#'
d_impact_block <-
  function(p_mean,
           p_var,
           etable,
           theta,
           mean_predictors,
           var_predictors,
           samp_size,
           num_items,
           num_quad,
           num_predictors) {

    # Obtain number of impact parameters.
    num_impact_parms <- length(p_mean) + length(p_var)

    # Make space for first and second derivatives.
    d1 <- matrix(0,nrow=num_impact_parms,ncol=1)
    d2 <- matrix(0,nrow=num_impact_parms,ncol=num_impact_parms)

    # Get latent mean and variance vectors.
    alpha <- mean_predictors %*% p_mean
    phi <- exp(var_predictors %*% p_var)

    # First and second derivatives for mean impact parameters.
    for(cov in 1:ncol(mean_predictors)) {
      d1_trace_mean <- vapply(1:num_quad,
                              function(x) {
                                mean_predictors[,cov]/phi*(theta[x]-alpha)
                              },numeric(samp_size))
      d2_trace_mean <- vapply(1:num_quad,
                              function(x) {
                                -mean_predictors[,cov]**2/phi
                              },numeric(samp_size))

      d1[cov,1] <- sum(etable*d1_trace_mean, na.rm = TRUE)
      d2[cov,cov] <- sum(etable*d2_trace_mean, na.rm = TRUE)
    }

    # First and second derivatives for variance impact parameters.
    for(cov in 1:ncol(var_predictors)) {
      d1_trace_var <-
        vapply(1:num_quad,
               function(x) {
                 .5*var_predictors[,cov]*(
                   (theta[x]-alpha)**2/phi - 1)
                 },numeric(samp_size))
      d2_trace_var <-
        vapply(1:num_quad,
               function(x) {
                 -.5*var_predictors[,cov]**2*(theta[x]-alpha)**2/phi
                 },numeric(samp_size))
      d1[ncol(mean_predictors)+cov,1] <-
        sum(etable*d1_trace_var, na.rm = TRUE)
      d2[ncol(mean_predictors)+cov,ncol(mean_predictors)+cov] <-
        sum(etable*d2_trace_var, na.rm = TRUE)
    }

    # Cross derivatives for mean and impact variance parameters.
    for(cov in 1:ncol(var_predictors)) {
      for(cov2 in 1:ncol(mean_predictors)) {

        d2_trace_cross <-
          vapply(1:num_quad,
                 function(x) {
                   var_predictors[,cov]*mean_predictors[,cov2]/phi*(
                     alpha-theta[x])
                 },numeric(samp_size))
        d2[ncol(mean_predictors)+cov,cov2] <-
          sum(etable*d2_trace_cross, na.rm = TRUE)

        if(cov2 > 1 && cov < cov2) {


          # Cross derivatives for mean parameters with different predictors.
          d2_trace_cross_mean <-
            vapply(1:num_quad,
                   function(x) {
                     -mean_predictors[,cov]*mean_predictors[,cov2]/phi
                   },numeric(samp_size))
          d2[cov2,cov] <-
            sum(etable*d2_trace_cross_mean,
                na.rm = TRUE)

          if(cov2 <= length(p_var)) {
          # Cross derivatives for variance parameters with different predictors.
          d2_trace_cross_var <-
            vapply(1:num_quad,
                   function(x) {
                     -.5*var_predictors[,cov]*var_predictors[,cov2]*(
                       theta[x]-alpha)**2/phi
                   },numeric(samp_size))
          d2[ncol(mean_predictors)+cov2,ncol(mean_predictors)+cov] <-
            sum(etable*d2_trace_cross_var,
                na.rm = TRUE)
          }
        }


      }
    }

    dlist <- list(d1,d2)

  }


#' Partial derivatives for mean and variance impact equation using observed score proxy.
#'
#' @param p_mean Vector of mean impact parameters.
#' @param p_var Vector of variance impact parameters.
#' @param prox_data Vector of observed proxy scores.
#' @param mean_predictors Possibly different matrix of predictors for the mean
#' impact equation.
#' @param var_predictors Possibly different matrix of predictors for the
#' variance impact equation.
#' @param samp_size Sample size in data set.
#' @param num_items Number of items in data set.
#' @param num_predictors Number of predictors in dataset.
#'
#' @return a \code{"list"} of first and second partial derivatives for impact equation (to
#' use with multivariate Newton-Rapshon and observed proxy scores)
#'
#' @keywords internal
#'
d_impact_block_proxy <-
  function(p_mean,
           p_var,
           prox_data,
           mean_predictors,
           var_predictors,
           samp_size,
           num_items,
           num_quad,
           num_predictors) {

    # Obtain number of impact parameters.
    num_impact_parms <- length(p_mean) + length(p_var)

    # Make space for first and second derivatives.
    d1 <- matrix(0,nrow=num_impact_parms,ncol=1)
    d2 <- matrix(0,nrow=num_impact_parms,ncol=num_impact_parms)

    # Get latent mean and variance vectors.
    alpha <- mean_predictors %*% p_mean
    phi <- exp(var_predictors %*% p_var)

    # First and second derivatives for mean impact parameters.
    for(cov in 1:ncol(mean_predictors)) {
      d1_trace_mean <- mean_predictors[,cov]/phi*(prox_data-alpha)
      d2_trace_mean <- -mean_predictors[,cov]**2/phi

      d1[cov,1] <- sum(d1_trace_mean, na.rm = TRUE)
      d2[cov,cov] <- sum(d2_trace_mean, na.rm = TRUE)
    }

    # First and second derivatives for variance impact parameters.
    for(cov in 1:ncol(var_predictors)) {
      d1_trace_var <- .5*var_predictors[,cov]*((prox_data-alpha)**2/phi - 1)
      d2_trace_var <- -.5*var_predictors[,cov]**2*(prox_data-alpha)**2/phi
      d1[ncol(mean_predictors)+cov,1] <- sum(d1_trace_var, na.rm = TRUE)
      d2[ncol(mean_predictors)+cov,ncol(mean_predictors)+cov] <- sum(d2_trace_var, na.rm = TRUE)
    }

    # Cross derivatives for mean and impact variance parameters.
    for(cov in 1:ncol(var_predictors)) {
      for(cov2 in 1:ncol(mean_predictors)) {

        d2_trace_cross <- var_predictors[,cov]*mean_predictors[,cov2]/phi*(alpha-prox_data)
        d2[ncol(mean_predictors)+cov,cov2] <- sum(d2_trace_cross, na.rm = TRUE)

        if(cov2 > 1 && cov < cov2) {


          # Cross derivatives for mean parameters with different predictors.
          d2_trace_cross_mean <- -mean_predictors[,cov]*mean_predictors[,cov2]/phi
          d2[cov2,cov] <- sum(d2_trace_cross_mean, na.rm = TRUE)

          if(cov2 <= length(p_var)) {
            # Cross derivatives for variance parameters with different predictors.
            d2_trace_cross_var <-
              -.5*var_predictors[,cov]*var_predictors[,cov2]*(prox_data-alpha)**2/phi
            d2[ncol(mean_predictors)+cov2,ncol(mean_predictors)+cov] <- sum(d2_trace_cross_var,
                                                                            na.rm = TRUE)
          }
        }


      }
    }

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
#' @return a \code{"list"} of first and second partial derivatives for Bernoulli item likelihood (to
#' use with coordinate descent and univariate Newton-Raphson)
#'
#' @keywords internal
#'
d_bernoulli <-
  function(parm,
           p_item,
           etable_item,
           theta,
           pred_data,
           cov,
           samp_size,
           num_items,
           num_quad) {

  if(parm == "c0"){
    eta_d <- matrix(1, nrow = samp_size, ncol = num_quad)
  } else if(parm == "a0"){
    eta_d <- t(matrix(theta,ncol=samp_size,nrow=num_quad))
  } else if(parm == "c1"){
    eta_d <- matrix(pred_data[,cov],
                    ncol = num_quad,
                    nrow = samp_size)
  } else if(parm == "a1"){
    eta_d <- matrix(pred_data[,cov],
                    ncol = num_quad,
                    nrow = samp_size)*t(matrix(theta,
                                               ncol=samp_size,
                                               nrow=num_quad))
  }

  traceline <- bernoulli_traceline_pts(p_item,
                                       theta,
                                       pred_data,
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

#' Partial derivatives for binary items by item-blocks.
#'
#' @param p_item Vector of item parameters.
#' @param etable E-table for item.
#' @param theta Matrix of adaptive theta values.
#' @param pred_data Matrix or dataframe of DIF and/or impact predictors.
#' @param item_data_current Vector of current item responses.
#' @param samp_size Sample size in dataset.
#' @param num_items Number of items in dataset.
#' @param num_predictors Number of predictors in dataset.
#' @param num_quad Number of quadrature points used for approximating the
#' latent variable.
#'
#' @return a \code{"list"} of first and second partial derivatives for Bernoulli item likelihood (to
#' use with multivariate Newton-Raphson)
#'
#' @keywords internal
#'
d_bernoulli_itemblock <-
  function(p_item,
           etable,
           theta,
           pred_data,
           item_data_current,
           samp_size,
           num_items,
           num_predictors,
           num_quad) {

    # Make space for first and second derivatives.
    d1 <- matrix(0,nrow=length(p_item),ncol=1)
    d2 <- matrix(0,nrow=length(p_item),ncol=length(p_item))

    # First derivative for linear predictor w.r.t. theta.
    eta_d_a0 <- t(matrix(theta,
                         ncol=samp_size,
                         nrow=num_quad))

    # Get item response function.
    traceline <- bernoulli_traceline_pts(p_item,
                                         theta,
                                         pred_data,
                                         samp_size)

    # Get posterior probabilities for each response.
    etable_item <- lapply(1:2, function(x) etable)

    # Obtain E-tables for each response category.
      for(resp in 1:2) {
        etable_item[[resp]][which(
          !(item_data_current == resp)), ] <- 0
      }

    # Calculate first and second base derivatives.
    d1_base <- traceline*(etable_item[[2]]/traceline -
                            etable_item[[2]] -
                            etable_item[[1]])
    d2_base <- (-traceline + traceline**2)*etable

    # First and second derivative for c0.
    d1[1,1] <- sum(d1_base, na.rm = TRUE) #d1
    d2[1,1] <- sum(d2_base, na.rm = TRUE) #d2

    # First and second derivative for a0.
    d1[2,1] <- sum(eta_d_a0*d1_base, na.rm = TRUE) #d1
    d2[2,2] <- sum(eta_d_a0**2*d2_base, na.rm = TRUE) #d2

    # Cross derivative for c0 and a0.
    d2[2,1] <- sum(eta_d_a0*d2_base, na.rm = TRUE) #d2

    # Cycle through predictors (outer cycle).
    for(cov in 1:num_predictors) {

      # First derivative for linear predictor w.r.t. covariate.
      cov_matrix <- matrix(pred_data[,cov],
                           ncol = num_quad,
                           nrow = samp_size)

      # First and second derivatives for c1.
      d1[2+cov,1] <-
        sum(cov_matrix*d1_base,
            na.rm = TRUE) #d1
      d2[2+cov,2+cov] <-
        sum(cov_matrix**2*d2_base,
            na.rm = TRUE) #d2

      # First and second derivatives for a1.
      d1[2+num_predictors+cov,1] <-
        sum(cov_matrix*eta_d_a0*d1_base,
            na.rm = TRUE) #d1
      d2[2+num_predictors+cov,2+num_predictors+cov] <-
        sum((cov_matrix*eta_d_a0)**2*d2_base,
            na.rm = TRUE) #d2

      # Cross derivatives for c0 and c1.
      d2[2+cov,1] <-
        sum(cov_matrix*d2_base,
            na.rm = TRUE) #d2

      # Cross derivatives for c0 and a1, as well as a0 and c1.
      d2[2+num_predictors+cov,1] <- d2[2+cov,2] <-
        sum(cov_matrix*eta_d_a0*d2_base,
            na.rm = TRUE) #d2

      # Cross derivatives for a0 and a1.
      d2[2+num_predictors+cov,2] <-
        sum(cov_matrix*eta_d_a0**2*d2_base,
            na.rm = TRUE) #d2

      # Cycle through predictors (inner cycle).
      for(cov2 in 1:num_predictors) {


        if(cov == cov2) {

          # Cross derivatives with same predictor for c1 and a1.
          d2[2+num_predictors+cov,2+cov2] <-
            sum(cov_matrix**2*eta_d_a0*d2_base,
                na.rm = TRUE) #d2

        } else {

          # First derivatives for linear predictor w.r.t. second covariate.
          cov2_matrix <- matrix(pred_data[,cov2],
                                ncol = num_quad,
                                nrow = samp_size)

          # Cross derivatives with different predictor for c1 and a1.
          d2[2+num_predictors+cov,2+cov2] <-
            sum(cov_matrix*cov2_matrix*eta_d_a0*d2_base,
                na.rm = TRUE) #d2

          if(cov2 > 1 && cov < cov2) {

            # Cross derivatives with different predictor for c1 and c1.
            d2[2+cov2,2+cov] <-
              sum(cov_matrix*cov2_matrix*d2_base, #d2
                  na.rm = TRUE)

            # Cross derivatives with different predictor for a1 and a1.
            d2[2+num_predictors+cov2,2+num_predictors+cov] <-
              sum(cov_matrix*cov2_matrix*eta_d_a0**2*d2_base, #a1a1
                  na.rm = TRUE)
          }
        }

      }

    }

    dlist <- list(d1,d2)

  }

#' Partial derivatives for binary items by item-blocks using observed score proxy.
#'
#' @param p_item Vector of item parameters.
#' @param pred_data Matrix or dataframe of DIF and/or impact predictors.
#' @param item_data_current Vector of current item responses.
#' @param prox_data Vector of observed proxy scores.
#' @param samp_size Sample size in dataset.
#' @param num_items Number of items in dataset.
#' @param num_predictors Number of predictors in dataset.
#'
#' @return a \code{"list"} of first and second partial derivatives for Bernoulli item likelihood (to
#' use with multivariate Newton-Raphson and observed proxy scores)
#'
#' @keywords internal
#'
d_bernoulli_itemblock_proxy <-
  function(p_item,
           pred_data,
           item_data_current,
           prox_data,
           samp_size,
           num_items,
           num_predictors,
           num_quad) {

    # Make space for first and second derivatives.
    d1 <- matrix(0,nrow=length(p_item),ncol=1)
    d2 <- matrix(0,nrow=length(p_item),ncol=length(p_item))

    # First derivative for linear predictor w.r.t. theta.
    eta_d_a0 <- prox_data

    # Get item response function.
    traceline <- bernoulli_traceline_pts_proxy(p_item,
                                               prox_data,
                                               pred_data)

    # Obtain tracelines for each response category.


    # Calculate first and second base derivatives.
    d1_base <- matrix(0, nrow = nrow(traceline), ncol = 1)
    d1_base[item_data_current == 1,] <- -traceline[item_data_current == 1,]
    d1_base[item_data_current == 2,] <- (1 - traceline)[item_data_current == 2,]

    d2_base <- -traceline + traceline**2

    # First and second derivative for c0.
    d1[1,1] <- sum(d1_base, na.rm = TRUE) #d1
    d2[1,1] <- sum(d2_base, na.rm = TRUE) #d2

    # First and second derivative for a0.
    d1[2,1] <- sum(eta_d_a0*d1_base, na.rm = TRUE) #d1
    d2[2,2] <- sum(eta_d_a0**2*d2_base, na.rm = TRUE) #d2

    # Cross derivative for c0 and a0.
    d2[2,1] <- sum(eta_d_a0*d2_base, na.rm = TRUE) #d2

    # Cycle through predictors (outer cycle).
    for(cov in 1:num_predictors) {

      # First derivative for linear predictor w.r.t. covariate.
      cov_matrix <- pred_data[,cov]

      # First and second derivatives for c1.
      d1[2+cov,1] <-
        sum(cov_matrix*d1_base,
            na.rm = TRUE) #d1
      d2[2+cov,2+cov] <-
        sum(cov_matrix**2*d2_base,
            na.rm = TRUE) #d2

      # First and second derivatives for a1.
      d1[2+num_predictors+cov,1] <-
        sum(cov_matrix*eta_d_a0*d1_base,
            na.rm = TRUE) #d1
      d2[2+num_predictors+cov,2+num_predictors+cov] <-
        sum((cov_matrix*eta_d_a0)**2*d2_base,
            na.rm = TRUE) #d2

      # Cross derivatives for c0 and c1.
      d2[2+cov,1] <-
        sum(cov_matrix*d2_base,
            na.rm = TRUE) #d2

      # Cross derivatives for c0 and a1, as well as a0 and c1.
      d2[2+num_predictors+cov,1] <- d2[2+cov,2] <-
        sum(cov_matrix*eta_d_a0*d2_base,
            na.rm = TRUE) #d2

      # Cross derivatives for a0 and a1.
      d2[2+num_predictors+cov,2] <-
        sum(cov_matrix*eta_d_a0**2*d2_base,
            na.rm = TRUE) #d2

      # Cycle through predictors (inner cycle).
      for(cov2 in 1:num_predictors) {


        if(cov == cov2) {

          # Cross derivatives with same predictor for c1 and a1.
          d2[2+num_predictors+cov,2+cov2] <-
            sum(cov_matrix**2*eta_d_a0*d2_base,
                na.rm = TRUE) #d2

        } else {

          # First derivatives for linear predictor w.r.t. second covariate.
          cov2_matrix <- pred_data[,cov2]

          # Cross derivatives with different predictor for c1 and a1.
          d2[2+num_predictors+cov,2+cov2] <-
            sum(cov_matrix*cov2_matrix*eta_d_a0*d2_base,
                na.rm = TRUE) #d2

          if(cov2 > 1 && cov < cov2) {

            # Cross derivatives with different predictor for c1 and c1.
            d2[2+cov2,2+cov] <-
              sum(cov_matrix*cov2_matrix*d2_base, #d2
                  na.rm = TRUE)

            # Cross derivatives with different predictor for a1 and a1.
            d2[2+num_predictors+cov2,2+num_predictors+cov] <-
              sum(cov_matrix*cov2_matrix*eta_d_a0**2*d2_base, #a1a1
                  na.rm = TRUE)
          }
        }

      }

    }

    dlist <- list(d1,d2)


  }

#' Partial derivatives for ordinal items.
#'
#' @param parm Item parameter being maximized.
#' @param p_item Vector of item parameters.
#' @param etable_item E-table for impact.
#' @param theta Matrix of adaptive theta values.
#' @param pred_data Matrix or dataframe of DIF and/or impact predictors.
#' @param thr Threshold value being maximized.
#' @param cov Covariate being maximized.
#' @param samp_size Sample size in dataset.
#' @param num_items Number of items in dataset.
#' @param num_quad Number of quadrature points used for approximating the
#' latent variable.
#'
#' @return a \code{"list"} of first and second partial derivatives for categorical item likelihood
#' (to use with coordinate descent and univariate Newton-Raphson)
#'
#' @keywords internal
#'
d_categorical <-
  function(parm,
           p_item,
           etable_item,
           theta,
           pred_data,
           thr,
           cov,
           samp_size,
           num_responses_item,
           num_items,
           num_quad) {

    if(parm == "c0"){
      eta_d <- matrix(1, nrow = samp_size, ncol = num_quad)
    } else if(parm == "a0"){
      eta_d <- t(matrix(theta,
                        ncol=samp_size,
                        nrow=num_quad))
    } else if(parm == "c1"){
      eta_d <- matrix(pred_data[,cov],
                      ncol = num_quad,
                      nrow = samp_size)
    } else if(parm == "a1"){
      eta_d <- matrix(pred_data[,cov],
                      ncol = num_quad,
                      nrow = samp_size)*t(matrix(theta,
                                                 ncol=samp_size,
                                                 nrow=num_quad))
    }

    cum_traceline <- cumulative_traceline_pts(p_item,
                                              theta,
                                              pred_data,
                                              samp_size,
                                              num_responses_item,
                                              num_quad)

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
      d2 <- sum(etable_item[[thr]] /
                  (cum_traceline[[thr-1]] -
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

#' Partial derivatives for ordinal items.
#'
#' @param parm Item parameter being maximized.
#' @param p_item Vector of item parameters.
#' @param etable E-table for impact.
#' @param theta Matrix of adaptive theta values.
#' @param pred_data Matrix or dataframe of DIF and/or impact predictors.
#' @param item_data_current Vector of current item responses.
#' @param samp_size Sample size in dataset.
#' @param num_items Number of items in dataset.
#' @param num_predictors Number of predictors in dataset.
#' @param num_quad Number of quadrature points used for approximating the
#' latent variable.
#'
#' @return a \code{"list"} of first and second partial derivatives for categorical item likelihood
#' (to use with multivariate Newton-Raphson)
#'
#' @keywords internal
#'
d_categorical_itemblock <-
  function(parm,
           p_item,
           etable,
           theta,
           pred_data,
           item_data_current,
           samp_size,
           num_responses_item,
           num_items,
           num_predictors,
           num_quad) {

    # Make space for first and second derivatives.
    d1 <- matrix(0,nrow=length(p_item),ncol=1)
    d2 <- matrix(0,nrow=length(p_item),ncol=length(p_item))

    # First derivative for linear predictor w.r.t. theta.
    eta_d_a0 <- t(matrix(theta,
                         ncol=samp_size,
                         nrow=num_quad))

    # Get item response function.
    cum_traceline <- cumulative_traceline_pts(p_item,
                                              theta,
                                              pred_data,
                                              samp_size,
                                              num_responses_item,
                                              num_quad)

    # Get posterior probabilities for each response.
    etable_item <- lapply(1:num_responses_item, function(x) etable)

    # Obtain E-tables for each response category.
    for(resp in 1:num_responses_item) {
      etable_item[[resp]][which(
        !(item_data_current == resp)), ] <- 0
    }


    # Calculate first and second base derivatives (non-thresholds).
    d1_base <- d2_base <- 0
    for(resp in 1:num_responses_item) {

      if(resp == 1) {
        d1_base <- d1_base +
          -etable_item[[1]]*cum_traceline[[1]]
        d2_base <- d2_base +
          -etable_item[[1]]*(cum_traceline[[1]]*(1-cum_traceline[[1]]))
      } else if(resp == num_responses_item) {
        d1_base <- d1_base +
          etable_item[[num_responses_item]]*(
            1-cum_traceline[[num_responses_item-1]])
        d2_base <- d2_base +
          etable_item[[num_responses_item]]*(
            -cum_traceline[[num_responses_item-1]]*(
              1-cum_traceline[[num_responses_item-1]]))
      } else {
        d1_base <- d1_base +
          etable_item[[resp]]*((1-cum_traceline[[resp]]) -
                                 cum_traceline[[resp-1]])
        d2_base <- d2_base +
          etable_item[[resp]]*(cum_traceline[[resp-1]]**2 +
                              cum_traceline[[resp]]**2 -
                              cum_traceline[[resp-1]] -
                              cum_traceline[[resp]])
      }

    }

    # First and second derivative for c0.
    d1[1,1] <- sum(d1_base, na.rm = TRUE) #d1
    d2[1,1] <- sum(d2_base, na.rm = TRUE) #d2

    # First and second derivative for a0.
    d1[num_responses_item,1] <-
      sum(eta_d_a0*d1_base, na.rm = TRUE) #d1
    d2[num_responses_item,num_responses_item] <-
      sum(eta_d_a0**2*d2_base, na.rm = TRUE) #d2

    # Cross derivative for c0 and a0.
    d2[num_responses_item,1] <- sum(eta_d_a0*d2_base, na.rm = TRUE) #d2

    # Cycle through predictors (outer cycle).
    for(cov in 1:num_predictors) {

      # First derivative for linear predictor w.r.t. covariate.
      cov_matrix <- matrix(pred_data[,cov],
                           ncol = num_quad,
                           nrow = samp_size)

      # First and second derivatives for c1.
      d1[num_responses_item+cov,1] <-
        sum(cov_matrix*d1_base,
            na.rm = TRUE) #d1
      d2[num_responses_item+cov,num_responses_item+cov] <-
        sum(cov_matrix**2*d2_base,
            na.rm = TRUE) #d2

      # First and second derivatives for a1.
      d1[num_responses_item+num_predictors+cov,1] <-
        sum(cov_matrix*eta_d_a0*d1_base,
            na.rm = TRUE) #d1
      d2[num_responses_item+num_predictors+cov,
         num_responses_item+num_predictors+cov] <-
        sum((cov_matrix*eta_d_a0)**2*d2_base,
            na.rm = TRUE) #d2

      # Cross derivatives for c0 and c1.
      d2[num_responses_item+cov,1] <-
        sum(cov_matrix*d2_base,
            na.rm = TRUE) #d2

      # Cross derivatives for c0 and a1, as well as a0 and c1.
      d2[num_responses_item+num_predictors+cov,1] <-
        d2[num_responses_item+cov,num_responses_item] <-
        sum(cov_matrix*eta_d_a0*d2_base,
            na.rm = TRUE) #d2

      # Cross derivatives for a0 and a1.
      d2[num_responses_item+num_predictors+cov,num_responses_item] <-
        sum(cov_matrix*eta_d_a0**2*d2_base,
            na.rm = TRUE) #d2

      # Cycle through predictors (inner cycle).
      for(cov2 in 1:num_predictors) {


        if(cov == cov2) {

          # Cross derivatives with same predictor for c1 and a1.
          d2[num_responses_item+num_predictors+cov,num_responses_item+cov2] <-
            sum(cov_matrix**2*eta_d_a0*d2_base,
                na.rm = TRUE) #d2

        } else {

          # First derivatives for linear predictor w.r.t. second covariate.
          cov2_matrix <- matrix(pred_data[,cov2],
                                ncol = num_quad,
                                nrow = samp_size)

          # Cross derivatives with different predictor for c1 and a1.
          d2[num_responses_item+num_predictors+cov,num_responses_item+cov2] <-
            sum(cov_matrix*cov2_matrix*eta_d_a0*d2_base,
                na.rm = TRUE) #d2

          if(cov2 > 1 && cov < cov2) {

            # Cross derivatives with different predictor for c1 and c1.
            d2[num_responses_item+cov2,num_responses_item+cov] <-
              sum(cov_matrix*cov2_matrix*d2_base, #d2
                  na.rm = TRUE)

            # Cross derivatives with different predictor for a1 and a1.
            d2[num_responses_item+num_predictors+cov2,
               num_responses_item+num_predictors+cov] <-
              sum(cov_matrix*cov2_matrix*eta_d_a0**2*d2_base, #a1a1
                  na.rm = TRUE)
          }
        }

      }

    }

    # Threshold derivatives.
    for(thr in 2:(num_responses_item-1)) {
      d1_base_thr <- cum_traceline[[thr]]*(1-cum_traceline[[thr]])
      d2_base_thr <- d1_base_thr*(1 - 2*cum_traceline[[thr]])

      cat_traceline1 <- cum_traceline[[thr-1]] - cum_traceline[[thr]]
      cat_traceline2 <- cum_traceline[[thr]]
      if(thr < (num_responses_item-1)) {
        cat_traceline2 <- cat_traceline2 - cum_traceline[[thr+1]]
      }

      d1[thr,1] <-
        sum(d1_base_thr*(-etable_item[[thr]] / cat_traceline1 +
                           etable_item[[thr+1]] / cat_traceline2),
            na.rm = TRUE)

      d2[thr,thr] <-
        sum(etable_item[[thr]] / cat_traceline1 *
              (d2_base_thr + d1_base_thr**2/cat_traceline1),
            na.rm = TRUE) -
        sum(etable_item[[thr+1]] / cat_traceline2 *
              (d2_base_thr - d1_base_thr**2/cat_traceline2),
            na.rm = TRUE)

      d2[thr,1] <-
        sum(etable_item[[thr]])


      # for(thr2 in 3:(num_responses_item-1))
      #
      # for(cov in 1:num_predictors) {
      #   # First derivative for linear predictor w.r.t. covariate.
      #   cov_matrix <- matrix(pred_data[,cov],
      #                        ncol = num_quad,
      #                        nrow = samp_size)
      # }

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
#' @param pred_data Matrix or dataframe of DIF and/or impact predictors.
#' @param cov Covariate being maximized.
#' @param samp_size Sample size in dataset.
#' @param num_items Number of items in dataset.
#' @param num_quad Number of quadrature points used for approximating the
#' latent variable.
#'
#' @return a \code{"list"} of first and second partial derivatives for mean value of Gaussian item
#' likelihood (to use with coordinate descent and univariate Newton-Raphson)
#'
#' @keywords internal
#'
d_mu_gaussian <-
  function(parm,
           p_item,
           etable_item,
           theta,
           responses_item,
           pred_data,
           cov,
           samp_size,
           num_items,
           num_quad) {

  if(parm == "c0"){
    eta_d <- matrix(1, nrow = samp_size, ncol = num_quad)
  } else if(parm == "a0"){
    eta_d <- t(replicate(n=samp_size, theta))
  } else if(parm == "c1"){
    eta_d <- matrix(rep(pred_data[,cov], num_quad),
                    ncol = num_quad,
                    nrow = samp_size)
  } else if(parm == "a1"){
    eta_d <- matrix(rep(pred_data[,cov], num_quad),
                    ncol = num_quad,
                    nrow = samp_size)*t(replicate(n=samp_size, theta))
  }


  # Get latent mean and variance vectors.
  mu <- sapply(theta,
              function(x) {
                (p_item[grep("c0",names(p_item),fixed=T)] +
                   pred_data %*%
                   p_item[grep("c1",names(p_item),fixed=T)]) +
                  (p_item[grep("a0",names(p_item),fixed=T)] +
                     pred_data %*%
                     p_item[grep("a1",names(p_item),fixed=T)])*x
                })
  sigma <- sqrt(p_item[grep("s0",names(p_item))][1]*exp(
    pred_data %*% p_item[grep("s1",names(p_item))]
    ))


  d1_trace <- t(sapply(1:samp_size,
                       function(x) {
                         eta_d[x,]/sigma[x]**2*(responses_item[x] - mu[x,])
                         }))
  d2_trace <- t(sapply(1:samp_size,
                       function(x) {
                         -eta_d[x,]**2 / sigma[x]**2
                         }))

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
#' @param pred_data Matrix or dataframe of DIF and/or impact predictors.
#' @param cov Covariate being maximized.
#' @param samp_size Sample size in dataset.
#' @param num_items Number of items in dataset.
#' @param num_quad Number of quadrature points used for approximating the
#' latent variable.
#'
#' @return a \code{"list"} of first and second partial derivatives for variance value of Gaussian
#' item likelihood (to use with coordinate descent and univariate Newton-Raphson)
#'
#' @keywords internal
#'
d_sigma_gaussian <-
  function(parm,
           p_item,
           etable_item,
           theta,
           responses_item,
           pred_data,
           cov,
           samp_size,
           num_items,
           num_quad) {

  sigma <- sqrt(p_item[grep("s0",names(p_item))][1]*exp(
    pred_data %*% p_item[grep("s1",names(p_item))]))
  mu <- sapply(theta,
              function(x) {
                (p_item[grep("c0",names(p_item),fixed=T)] +
                   pred_data %*%
                   p_item[grep("c1",names(p_item),fixed=T)]) +
                  (p_item[grep("a0",names(p_item),fixed=T)] +
                     pred_data %*%
                     p_item[grep("a1",names(p_item),fixed=T)])*x
                })

  if(parm == "s0") {
    eta_d1 <- sapply(1:samp_size,
                     function(x) {
                       exp(pred_data[x,] %*%
                             p_item[grep("s1",names(p_item))]) / (2*sigma[x])
                       })
    eta_d2 <- sapply(1:samp_size,
                     function(x) {
                       -exp(pred_data[x,] %*%
                              p_item[grep("s1",names(p_item))])**2 /
                         (4*sigma[x]**3)
                       })
  } else if(parm == "s1") {
    eta_d1 <- sapply(1:samp_size,
                     function(x) {
                       sigma[x]*pred_data[x,cov] / 2
                       })
    eta_d2 <- sapply(1:samp_size,
                     function(x) {
                       sigma[x]*pred_data[x,cov]**2 / 4
                       })
  }


  d1_trace <- t(sapply(1:samp_size,
                       function(x) {
                         eta_d1[x]*((responses_item[x]-mu[x,])**2 /
                                      sigma[x]**3 -
                                      1/sigma[x])
                         }))

  if(parm == "s0") {
    d2_trace <- t(sapply(1:samp_size,
                         function(x) {
                           eta_d1[x]**2*(1 / sigma[x]**2 -
                                           3*(responses_item[x] - mu[x,])**2 /
                                           sigma[x]**4) +
                             eta_d2[x]*((responses_item[x] - mu[x,])**2 /
                                          sigma[x]**3 - 1/sigma[x])
                           }))
  } else if(parm == "s1") {
    d2_trace <- t(sapply(1:samp_size,
                         function(x) {
                           -2*eta_d2[x]*(sigma[x]**(-3)*(responses_item[x] -
                                                           mu[x,])**2)
                           }))
  }

  d1 <- sum(etable_item[[1]]*d1_trace, na.rm = TRUE)
  d2 <- sum(etable_item[[1]]*d2_trace, na.rm = TRUE)

  dlist <- list(d1,d2)

}



