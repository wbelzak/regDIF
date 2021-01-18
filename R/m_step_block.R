#' Maximization step using latent variable and item response blocks.
#'
#' @param p List of parameters.
#' @param item_data Matrix or data frame of item responses.
#' @param pred_data Matrix or data frame of DIF and/or impact predictors.
#' @param mean_predictors Possibly different matrix of predictors for the mean
#' impact equation.
#' @param var_predictors Possibly different matrix of predictors for the
#' variance impact equation.
#' @param etable E-table matrix for item and impact equations, in addition to
#' theta values (possibly adaptive).
#' @param item_type Optional character value or vector indicating the type of
#' item to be modeled.
#' @param pen_type Character value indicating the penalty function to use.
#' @param tau_current A single numeric value of tau that exists within
#' \code{tau_vec}.
#' @param alpha Numeric value indicating the alpha parameter in the elastic net
#' penalty function.
#' @param gamma Numeric value indicating the gamma parameter in the MCP
#' function.
#' @param anchor Optional numeric value or vector indicating which item
#' response(s) are anchors (e.g., \code{anchor = 1}).
#' @param samp_size Sample size in data set.
#' @param num_responses Number of responses for each item.
#' @param num_items Number of items in data set.
#' @param num_quad Number of quadrature points used for approximating the
#' latent variable.
#' @param num_predictors Number of predictors.
#'
#' @keywords internal
#'
Mstep_block <-
  function(p,
           item_data,
           pred_data,
           mean_predictors,
           var_predictors,
           etable,
           item_type,
           pen_type,
           tau_current,
           alpha,
           gamma,
           anchor,
           samp_size,
           num_responses,
           num_items,
           num_quad,
           num_predictors) {

    # Latent impact updates.
    anl_deriv_impact <- d_impact_block(p[[num_items+1]],
                                       p[[num_items+2]],
                                       etable$etable,
                                       etable$theta,
                                       mean_predictors,
                                       var_predictors,
                                       samp_size,
                                       num_items,
                                       num_quad,
                                       num_predictors)
    m <- c(p[[num_items+1]],p[[num_items+2]]) -
      solve(anl_deriv_impact[[2]]) %*% anl_deriv_impact[[1]]

    for(cov in 1:ncol(mean_predictors)){
      p[[num_items+1]][[cov]] <- m[cov]
    }
    for(cov in 1:ncol(var_predictors)){
      p[[num_items+2]][[cov]] <- m[ncol(mean_predictors)+cov]
    }


    # Item response updates.
    for (item in 1:num_items) {
    # foreach(item=1:num_items) %do% {

      # Bernoulli responses.
      if(num_responses[item] == 2) {

        anl_deriv_item <- d_bernoulli_itemblock(p[[item]],
                                                etable$etable,
                                                etable$theta,
                                                pred_data,
                                                item_data[,item],
                                                samp_size,
                                                num_items,
                                                num_predictors,
                                                num_quad)

        z <- p[[item]] - solve(anl_deriv_item[[2]]) %*% anl_deriv_item[[1]]

        # c0 update.
        p[[item]][[1]] <- z[1]

        # a0 update.
        if(item_type[item] != "Rasch") p[[item]][[2]] <- z[2]

        # Don't update DIF estimates if anchor item.
        if(any(item == anchor)) next

        p2 <- unlist(p)

        for(cov in 1:num_predictors) {

          # End routine if only one anchor item is left on each covariate
          # for each item parameter.
          if(is.null(anchor) &&
             sum(p2[grep(paste0("c1(.*?)cov",cov),names(p2))] != 0) >
             (num_items - 2) &&
             alpha == 1){
            next
          }

          # c1 updates.
          p[[item]][[2+cov]] <-
            ifelse(pen_type == "lasso",
                   soft_threshold(z[2+cov],
                                   alpha,
                                   tau_current),
                   firm_threshold(z[2+cov],
                                   alpha,
                                   tau_current,
                                   gamma))

          # a1 updates.
          if(item_type[item] == "Rasch") next

          # End routine if only one anchor item is left on each covariate
          # for each item parameter.
          if(is.null(anchor) &&
             sum(p2[grep(paste0("a1(.*?)cov",cov),names(p2))] != 0) >
             (num_items - 2) &&
             alpha == 1){
            next
          }

          p[[item]][[2+num_predictors+cov]] <-
            ifelse(pen_type == "lasso",
                   soft_threshold(z[2+num_predictors+cov],
                                   alpha,
                                   tau_current),
                   firm_threshold(z[2+num_predictors+cov],
                                   alpha,
                                   tau_current,
                                   gamma))

        }


      } else if(num_responses[item] > 2) {

        # Get posterior probabilities for each response.
        etable_item <- replicate(n=num_responses[item],
                                 etable$etable,
                                 simplify = F)

        # Obtain E-tables for each response category.
        for(resp in 1:num_responses[item]) {
          etable_item[[resp]][which(
            !(item_data[,item] == resp)), ] <- 0
        }

        # Intercept updates.
        anl_deriv <- d_categorical("c0",
                                   p[[item]],
                                   etable_item,
                                   etable$theta,
                                   pred_data,
                                   thr=-1,
                                   cov=-1,
                                   samp_size,
                                   num_responses[[item]],
                                   num_items,
                                   num_quad)
        p[[item]][[1]] <- p[[item]][[1]] - anl_deriv[[1]]/anl_deriv[[2]]

        # Threshold updates.
        for(thr in 2:(num_responses[item]-1)) {
          anl_deriv <- d_categorical("c0",
                                     p[[item]],
                                     etable_item,
                                     etable$theta,
                                     pred_data,
                                     thr=thr,
                                     cov=-1,
                                     samp_size,
                                     num_responses[[item]],
                                     num_items,
                                     num_quad)
          p[[item]][[thr]] <- p[[item]][[thr]] - anl_deriv[[1]]/anl_deriv[[2]]
        }

        # Slope updates.
        if(item_type[item] != "Rasch") {
          anl_deriv <- d_categorical("a0",
                                     p[[item]],
                                     etable_item,
                                     etable$theta,
                                     pred_data,
                                     thr=-1,
                                     cov=-1,
                                     samp_size,
                                     num_responses[[item]],
                                     num_items,
                                     num_quad)
          p[[item]][[num_responses[[item]]]] <-
            p[[item]][[num_responses[[item]]]] - anl_deriv[[1]]/anl_deriv[[2]]
        }

        if(any(item == anchor)) next

        p2 <- unlist(p)

        # Intercept DIF updates.
        for(cov in 1:num_predictors) {

          # End routine if only one anchor item is left on each covariate
          # for each item parameter.
          if(is.null(anchor) &
             sum(p2[grep(paste0("c1(.*?)cov",cov),names(p2))] != 0) >
             (num_items - 2) &&
             alpha == 1){
            next
          }

          anl_deriv <- d_categorical("c1",
                                     p[[item]],
                                     etable_item,
                                     etable$theta,
                                     pred_data,
                                     thr=-1,
                                     cov,
                                     samp_size,
                                     num_responses[[item]],
                                     num_items,
                                     num_quad)
          z <- p[[item]][[num_responses[[item]]+cov]] -
            anl_deriv[[1]]/anl_deriv[[2]]
          p[[item]][[num_responses[[item]]+cov]] <-
            ifelse(pen_type == "lasso",
                   soft_thresh_cpp(z,alpha,tau_current),
                   firm_thresh_cpp(z,alpha,tau_current,gamma))
        }

        # Slope DIF updates.
        for(cov in 1:num_predictors) {

          # End routine if only one anchor item is left on each covariate
          # for each item parameter.
          if(is.null(anchor) &
             sum(p2[grep(paste0("a1(.*?)cov",cov),names(p2))] != 0) >
             (num_items - 2) &
             alpha == 1){
            next
          }

          if(item_type[item] != "Rasch") {
            anl_deriv <- d_categorical("a1",
                                       p[[item]],
                                       etable_item,
                                       etable$theta,
                                       pred_data,
                                       thr=-1,
                                       cov,
                                       samp_size,
                                       num_responses[[item]],
                                       num_items,
                                       num_quad)
            z <- p[[item]][[length(p[[item]])-ncol(pred_data)+cov]] -
              anl_deriv[[1]]/anl_deriv[[2]]
            p[[item]][[length(p[[item]])-ncol(pred_data)+cov]] <-
              ifelse(pen_type == "lasso",
                     soft_thresh_cpp(z,alpha,tau_current),
                     firm_thresh_cpp(z,alpha,tau_current,gamma))
          }
        }



        # Gaussian responses.
      } else if(num_responses[item] == 1) {

        # Intercept updates.
        anl_deriv <- d_mu_gaussian("c0",
                                   p[[item]],
                                   etable$etable,
                                   etable$theta,
                                   item_data[,item],
                                   pred_data,
                                   cov=NULL,
                                   samp_size,
                                   num_items,
                                   num_quad)
        p[[item]][[1]] <- p[[item]][[1]] - anl_deriv[[1]]/anl_deriv[[2]]

        # Slope updates.
        if(item_type[item] != "Rasch") {
          a0_parms <- grep(paste0("a0_itm",item,"_"),names(p[[item]]),fixed=T)
          anl_deriv <- d_mu_gaussian("a0",
                                     p[[item]],
                                     etable$etable,
                                     etable$theta,
                                     item_data[,item],
                                     pred_data,
                                     cov=NULL,
                                     samp_size,
                                     num_items,
                                     num_quad)
          p[[item]] <- p[[item]][a0_parms] - anl_deriv[[1]]/anl_deriv[[2]]
        }

        # Residual updates.
        s0_parms <- grep(paste0("s0_itm",item,"_"),names(p[[item]]),fixed=T)
        anl_deriv <- d_sigma_gaussian("s0",
                                      p[[item]],
                                      etable$etable,
                                      etable$theta,
                                      item_data[,item],
                                      pred_data,
                                      cov=NULL,
                                      samp_size,
                                      num_items,
                                      num_quad)
        p[[item]][s0_parms][1] <- p[[item]][s0_parms][1] -
          anl_deriv[[1]]/anl_deriv[[2]]


        if(!any(item == anchor)) {

          # Residual DIF updates.
          for(cov in 1:num_predictors) {
            s1_parms <-
              grep(paste0("s1_itm",item,"_cov",cov),names(p[[item]]),fixed=T)
            anl_deriv <- d_sigma_gaussian("s1",
                                          p[[item]],
                                          etable$etable,
                                          etable$theta,
                                          item_data[,item],
                                          pred_data,
                                          cov=cov,
                                          samp_size,
                                          num_items,
                                          num_quad)
            p[[item]][s1_parms][1] <- p[[item]][s1_parms][1] -
              anl_deriv[[1]]/anl_deriv[[2]]
          }

          p2 <- unlist(p)

          # Intercept DIF updates.
          for(cov in 1:num_predictors){

            # End routine if only one anchor item is left on each covariate
            # for each item parameter.
            if(is.null(anchor) &
               sum(p2[grep(paste0("c1(.*?)cov",cov),names(p2))] != 0) >
               (num_items - 2) &
               alpha == 1){
              next
            }

            c1_parms <-
              grep(paste0("c1_itm",item,"_cov",cov),names(p[[item]]),fixed=T)
            anl_deriv <- d_mu_gaussian("c1",
                                       p[[item]],
                                       etable$etable,
                                       etable$theta,
                                       item_data[,item],
                                       pred_data,
                                       cov,
                                       samp_size,
                                       num_items,
                                       num_quad)
            z <- p[[item]][c1_parms] - anl_deriv[[1]]/anl_deriv[[2]]
            p[[item]][c1_parms] <-
              ifelse(pen_type == "lasso",
                     soft_thresh_cpp(z,alpha,tau_current),
                     firm_thresh_cpp(z,alpha,tau_current,gamma))
          }

          # Slope DIF updates.
          for(cov in 1:num_predictors){

            # End routine if only one anchor item is left on each covariate
            # for each item parameter.
            if(is.null(anchor) &
               sum(p2[grep(paste0("a1(.*?)cov",cov),names(p2))] != 0) >
               (num_items - 2) &
               alpha == 1){
              next
            }

            if(item_type[item] != "Rasch"){
              a1_parms <-
                grep(paste0("a1_itm",item,"_cov",cov),names(p[[item]]),fixed=T)
              anl_deriv <- d_mu_gaussian("a1",
                                         p[[item]],
                                         etable$etable,
                                         etable$theta,
                                         item_data[,item],
                                         pred_data,
                                         cov,
                                         samp_size,
                                         num_items,
                                         num_quad)
              z <- p[[item]][a1_parms] - anl_deriv[[1]]/anl_deriv[[2]]
              p[[item]][a1_parms] <- ifelse(pen_type == "lasso",
                                            soft_thresh_cpp(z,alpha,tau_current),
                                            firm_thresh_cpp(z,alpha,tau_current,gamma))
            }
          }
        }

      }

    }

    return(p)

  }

#' Final M-step using block estimation for identifying the minimum value of tau.
#'
#' @param p List of parameters.
#' @param item_data Matrix or data frame of item responses.
#' @param pred_data Matrix or data frame of DIF and/or impact predictors.
#' @param mean_predictors Possibly different matrix of predictors for the mean
#' impact equation.
#' @param var_predictors Possibly different matrix of predictors for the
#' variance impact equation.
#' @param etable Etable for item and impact equations, in addition to
#' theta values.
#' @param item_type Optional character value or vector indicating the type of
#' item to be modeled.
#' @param pen_type Character value indicating the penalty function to use.
#' @param tau_current A single numeric value of tau that exists within
#' \code{tau_vec}.
#' @param alpha Numeric value indicating the alpha parameter in the elastic net
#' penalty function.
#' @param gamma Numeric value indicating the gamma parameter in the MCP
#' function.
#' @param anchor Optional numeric value or vector indicating which item
#' response(s) are anchors (e.g., \code{anchor = 1}).
#' @param samp_size Sample size in data set.
#' @param num_responses Number of responses for each item.
#' @param num_items Number of items in data set.
#' @param num_quad Number of quadrature points used for approximating the
#' latent variable.
#' @param num_predictors Number of predictors.
#'
#' @keywords internal
#'
Mstep_block_idtau <-
  function(p,
           item_data,
           pred_data,
           mean_predictors,
           var_predictors,
           etable,
           item_type,
           pen_type,
           tau_current,
           alpha,
           gamma,
           anchor,
           samp_size,
           num_responses,
           num_items,
           num_quad,
           num_predictors) {


    # Last Mstep
    id_max_z <- 0

    # Update theta and etable.
    theta <- etable$theta
    etable <- etable$etable

    # Maximize independent logistic regressions.
    for (item in 1:num_items) {



      # Bernoulli responses.
      if(num_responses[item] == 2) {

        anl_deriv_item <- d_bernoulli_itemblock(p[[item]],
                                                etable,
                                                theta,
                                                pred_data,
                                                item_data[,item],
                                                samp_size,
                                                num_items,
                                                num_predictors,
                                                num_quad)

        z <- p[[item]] - solve(anl_deriv_item[[2]]) %*% anl_deriv_item[[1]]

        if(any(item == anchor)) next

        p2 <- unlist(p)

        # Intercept DIF updates.
        for(cov in 1:num_predictors) {

          # End routine if only one anchor item is left on each covariate
          # for each item parameter.
          if(is.null(anchor) &&
             sum(p2[grep(paste0("c1(.*?)cov",cov),names(p2))] != 0) >
             (num_items - 2) &&
             sum(p2[grep(paste0("a1(.*?)cov",cov),names(p2))] != 0) >
             (num_items - 2) &&
             alpha == 1){
            next
          }

          id_max_z <- c(id_max_z,
                        z[2+cov],
                        z[2+num_predictors+cov])
        }




      } else if(num_responses[item] > 2) {

        # Get posterior probabilities.
        etable_item <- replicate(n=num_responses[item],
                                 etable,
                                 simplify = F)

        # Obtain E-tables for each response category.
        if(num_responses[item] > 1) {
          for(resp in 1:num_responses[item]) {
            etable_item[[resp]][which(
              !(item_data[,item] == resp)), ] <- 0
          }
        }

        if(!any(item == anchor)){

          p2 <- unlist(p)

          # Intercept DIF updates.
          for(cov in 1:num_predictors) {

            # End routine if only one anchor item is left on each covariate
            # for each item parameter.
            if(is.null(anchor) &
               sum(p2[grep(paste0("c1(.*?)cov",cov),names(p2))] != 0) >
               (num_items - 1) &&
               alpha == 1){
              next
            }

            anl_deriv <- d_categorical("c1",
                                       p[[item]],
                                       etable_item,
                                       theta,
                                       pred_data,
                                       thr=-1,
                                       cov,
                                       samp_size,
                                       num_responses[[item]],
                                       num_items,
                                       num_quad)
            z <- p[[item]][[num_responses[[item]]+cov]] -
              anl_deriv[[1]]/anl_deriv[[2]]
            id_max_z <- c(id_max_z,z)
          }

          # Slope DIF updates.
          for(cov in 1:num_predictors) {

            # End routine if only one anchor item is left on each covariate
            # for each item parameter.
            if(is.null(anchor) &
               sum(p2[grep(paste0("a1(.*?)cov",cov),names(p2))] != 0) >
               (num_items - 1) &
               alpha == 1){
              next
            }

            if(item_type[item] != "Rasch") {
              anl_deriv <- d_categorical("a1",
                                         p[[item]],
                                         etable_item,
                                         theta,
                                         pred_data,
                                         thr=-1,
                                         cov,
                                         samp_size,
                                         num_responses[[item]],
                                         num_items,
                                         num_quad)
              z <- p[[item]][[length(p[[item]])-ncol(pred_data)+cov]] -
                anl_deriv[[1]]/anl_deriv[[2]]
              id_max_z <- c(id_max_z,z)
            }
          }

        }


        # Gaussian responses.
      } else if(num_responses[item] == 1) {

        if(!any(item == anchor)) {

          p2 <- unlist(p)

          # Intercept DIF updates.
          for(cov in 1:num_predictors){

            # End routine if only one anchor item is left on each covariate
            # for each item parameter.
            if(is.null(anchor) &
               sum(p2[grep(paste0("c1(.*?)cov",cov),names(p2))] != 0) >
               (num_items - 2) &
               alpha == 1){
              next
            }

            c1_parms <-
              grep(paste0("c1_itm",item,"_cov",cov),names(p[[item]]),fixed=T)
            anl_deriv <- d_mu_gaussian("c1",
                                       p[[item]],
                                       etable_item,
                                       theta,
                                       item_data[,item],
                                       pred_data,
                                       cov,
                                       samp_size,
                                       num_items,
                                       num_quad)
            z <- p[[item]][[c1_parms]] - anl_deriv[[1]]/anl_deriv[[2]]
            id_max_z <- c(id_max_z,z)
          }

          # Slope DIF updates.
          for(cov in 1:num_predictors){

            # End routine if only one anchor item is left on each covariate
            # for each item parameter.
            if(is.null(anchor) &
               sum(p2[grep(paste0("a1(.*?)cov",cov),names(p2))] != 0) >
               (num_items - 2) &
               alpha == 1){
              next
            }

            if(item_type[item] != "Rasch"){
              a1_parms <-
                grep(paste0("a1_itm",item,"_cov",cov),names(p[[item]]),fixed=T)
              anl_deriv <- d_mu_gaussian("a1",
                                         p[[item]],
                                         etable_item,
                                         theta,
                                         item_data[,item],
                                         pred_data,
                                         cov,
                                         samp_size,
                                         num_items,
                                         num_quad)
              z <- p[[item]][[a1_parms]] - anl_deriv[[1]]/anl_deriv[[2]]
              id_max_z <- c(id_max_z,z)
            }
          }
        }

      }



    }


    id_max_z <- max(abs(id_max_z))
    return(id_max_z)

  }
