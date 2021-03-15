#' Maximization step using coordinate descent optimization.
#'
#' @param p List of parameters.
#' @param item_data Matrix or data frame of item responses.
#' @param pred_data Matrix or data frame of DIF and/or impact predictors.
#' @param mean_predictors Possibly different matrix of predictors for the mean
#' impact equation.
#' @param var_predictors Possibly different matrix of predictors for the
#' variance impact equation.
#' @param eout E step output, including matrix for item and impact equations,
#' in addition to theta values (possibly adaptive).
#' @param item_type Optional character value or vector indicating the type of
#' item to be modeled.
#' @param pen_type Character value indicating the penalty function to use.
#' @param tau_current A single numeric value of tau that exists within
#' \code{tau_vec}.
#' @param pen Current penalty index.
#' @param alpha Numeric value indicating the alpha parameter in the elastic net
#' penalty function.
#' @param gamma Numeric value indicating the gamma parameter in the MCP
#' function.
#' @param anchor Optional numeric value or vector indicating which item
#' response(s) are anchors (e.g., \code{anchor = 1}).
#' @param final_control Control parameters.
#' @param samp_size Sample size in data set.
#' @param num_responses Number of responses for each item.
#' @param num_items Number of items in data set.
#' @param num_quad Number of quadrature points used for approximating the
#' latent variable.
#' @param num_predictors Number of predictors.
#' @param num_tau Logical indicating whether the minimum tau value needs to be
#' identified during the regDIF procedure.
#' @param max_tau Logical indicating whether to output the maximum tau value
#' needed to remove all DIF from the model.
#'
#' @keywords internal
#'
Mstep_cd2 <-
  function(p,
           item_data,
           pred_data,
           mean_predictors,
           var_predictors,
           eout,
           item_type,
           pen_type,
           tau_current,
           pen,
           alpha,
           gamma,
           anchor,
           final_control,
           samp_size,
           num_responses,
           num_items,
           num_quad,
           num_predictors,
           num_tau,
           max_tau) {

    # Set under-identified model to FALSE until proven TRUE.
    under_identified <- FALSE

    # Update theta and etable.
    theta <- eout$theta
    etable <- eout$etable
    p_cd <- p

    # Last Mstep
    if(max_tau) id_max_z <- 0


    # CD Maximization and print settings.
    lastp_cd_all <- p_cd
    eps_cd_all <- Inf
    iter_cd_all <- 1

    # Loop until convergence or maximum number of iterations.
    while(eps_cd_all > final_control$tol){



    # Latent mean impact updates.
    for(cov in 1:ncol(mean_predictors)) {

      # CD Maximization and print settings.
      lastp_cd <- p_cd
      eps_cd <- Inf
      iter_cd <- 1

      # Loop until convergence or maximum number of iterations.
      while(eps_cd > final_control$tol){

      anl_deriv <- d_alpha(c(p_cd[[num_items+1]],p_cd[[num_items+2]]),
                           etable,
                           theta,
                           mean_predictors,
                           var_predictors,
                           cov=cov,
                           samp_size,
                           num_items,
                           num_quad)
      p_cd[[num_items+1]][[cov]] <-
        p_cd[[num_items+1]][[cov]] - anl_deriv[[1]]/anl_deriv[[2]]

      # Update and check for convergence: Calculate the difference
      # in parameter estimates from current to previous.
      eps_cd = sqrt(sum((unlist(p_cd)-unlist(lastp_cd))^2))

      # Update parameter list.
      lastp_cd <- p_cd

      # Update the iteration number.
      iter_cd = iter_cd + 1


      }


    }

    # Latent variance impact updates.
    for(cov in 1:ncol(var_predictors)) {

      # CD Maximization and print settings.
      lastp_cd <- p_cd
      eps_cd <- Inf
      iter_cd <- 1

      # Loop until convergence or maximum number of iterations.
      while(eps_cd > final_control$tol){

      anl_deriv <- d_phi(c(p_cd[[num_items+1]],p_cd[[num_items+2]]),
                         etable,
                         theta,
                         mean_predictors,
                         var_predictors,
                         cov=cov,
                         samp_size,
                         num_items,
                         num_quad)
      p_cd[[num_items+2]][[cov]] <-
        p_cd[[num_items+2]][[cov]] - anl_deriv[[1]]/anl_deriv[[2]]

      # Update and check for convergence: Calculate the difference
      # in parameter estimates from current to previous.
      eps_cd = sqrt(sum((unlist(p_cd)-unlist(lastp_cd))^2))

      # Update parameter list.
      lastp_cd <- p_cd

      # Update the iteration number.
      iter_cd = iter_cd + 1
      }

    }


    # Item response updates.
    for (item in 1:num_items) {

      # Get posterior probabilities.
      etable_item <- lapply(1:num_responses[item], function(x) etable)

      # Obtain E-tables for each response category.
      if(num_responses[item] > 1) {
        for(resp in 1:num_responses[item]) {
          etable_item[[resp]][which(
            !(item_data[,item] == resp)), ] <- 0
        }
      }

      # Bernoulli responses.
      if(num_responses[item] == 2) {

        # CD Maximization and print settings.
        lastp_cd <- p_cd
        eps_cd <- Inf
        iter_cd <- 1

        # Loop until convergence or maximum number of iterations.
        while(eps_cd > final_control$tol){

        # Intercept updates.
        anl_deriv <- d_bernoulli("c0",
                                 p_cd[[item]],
                                 etable_item,
                                 theta,
                                 pred_data,
                                 cov=0,
                                 samp_size,
                                 num_items,
                                 num_quad)
        p_cd[[item]][[1]] <- p_cd[[item]][[1]] - anl_deriv[[1]]/anl_deriv[[2]]

        # Update and check for convergence: Calculate the difference
        # in parameter estimates from current to previous.
        eps_cd = sqrt(sum((unlist(p_cd)-unlist(lastp_cd))^2))

        # Update parameter list.
        lastp_cd <- p_cd

        # Update the iteration number.
        iter_cd = iter_cd + 1

        }

        # Slope updates.
        if(item_type[item] != "Rasch") {

          # CD Maximization and print settings.
          lastp_cd <- p_cd
          eps_cd <- Inf
          iter_cd <- 1

          # Loop until convergence or maximum number of iterations.
          while(eps_cd > final_control$tol){
          anl_deriv <- d_bernoulli("a0",
                                   p_cd[[item]],
                                   etable_item,
                                   theta,
                                   pred_data,
                                   cov=0,
                                   samp_size,
                                   num_items,
                                   num_quad)
          p_cd[[item]][[2]] <- p_cd[[item]][[2]] - anl_deriv[[1]]/anl_deriv[[2]]
          # Update and check for convergence: Calculate the difference
          # in parameter estimates from current to previous.
          eps_cd = sqrt(sum((unlist(p_cd)-unlist(lastp_cd))^2))

          # Update parameter list.
          lastp_cd <- p_cd

          # Update the iteration number.
          iter_cd = iter_cd + 1
          }

        }


        if(!any(item == anchor)) {

          p2_cd <- unlist(p_cd)

          # Intercept DIF updates.
          for(cov in 1:num_predictors) {

            # End routine if only one anchor item is left on each covariate
            # for each item parameter.
            if(is.null(anchor) &
               sum(p2_cd[grep(paste0("c1(.*?)cov",cov),names(p2_cd))] != 0) >
               (num_items - 1) &
               alpha == 1 &&
               (length(final_control$start.values) == 0 || pen > 1) &&
               num_tau >= 10){
              under_identified <- TRUE
              break
            }

            # CD Maximization and print settings.
            lastp_cd <- p_cd
            eps_cd <- Inf
            iter_cd <- 1

            # Loop until convergence or maximum number of iterations.
            while(eps_cd > final_control$tol){

            anl_deriv <- d_bernoulli("c1",
                                     p_cd[[item]],
                                     etable_item,
                                     theta,
                                     pred_data,
                                     cov,
                                     samp_size,
                                     num_items,
                                     num_quad)
            z <- p_cd[[item]][[2+cov]] - anl_deriv[[1]]/anl_deriv[[2]]
            if(max_tau) id_max_z <- c(id_max_z,z)
            p_cd[[item]][[2+cov]] <- ifelse(pen_type == "lasso",
                                         soft_thresh_cpp(z,alpha,tau_current),
                                         firm_thresh_cpp(z,alpha,tau_current,gamma))

            # Update and check for convergence: Calculate the difference
            # in parameter estimates from current to previous.
            eps_cd = sqrt(sum((unlist(p_cd)-unlist(lastp_cd))^2))

            # Update parameter list.
            lastp_cd <- p_cd

            # Update the iteration number.
            iter_cd = iter_cd + 1
            }

          }

          # Slope DIF updates.
          for(cov in 1:num_predictors) {

            # End routine if only one anchor item is left on each covariate
            # for each item parameter.
            if(is.null(anchor) &
               sum(p2_cd[grep(paste0("a1(.*?)cov",cov),names(p2_cd))] != 0) >
               (num_items - 1) &
               alpha == 1 &&
               (length(final_control$start.values) == 0 || pen > 1) &&
               num_tau >= 10){
              under_identified <- TRUE
              break
            }

            if(item_type[item] != "Rasch") {

              # CD Maximization and print settings.
              lastp_cd <- p_cd
              eps_cd <- Inf
              iter_cd <- 1

              # Loop until convergence or maximum number of iterations.
              while(eps_cd > final_control$tol){

              anl_deriv <- d_bernoulli("a1",
                                       p_cd[[item]],
                                       etable_item,
                                       theta,
                                       pred_data,
                                       cov,
                                       samp_size,
                                       num_items,
                                       num_quad)
              z <- p_cd[[item]][[2+num_predictors+cov]] -
                anl_deriv[[1]]/anl_deriv[[2]]
              if(max_tau) id_max_z <- c(id_max_z,z)
              p_cd[[item]][[2+num_predictors+cov]] <-
                ifelse(pen_type == "lasso",
                       soft_thresh_cpp(z,alpha,tau_current),
                       firm_thresh_cpp(z,alpha,tau_current,gamma))

              # Update and check for convergence: Calculate the difference
              # in parameter estimates from current to previous.
              eps_cd = sqrt(sum((unlist(p_cd)-unlist(lastp_cd))^2))

              # Update parameter list.
              lastp_cd <- p_cd

              # Update the iteration number.
              iter_cd = iter_cd + 1
              }

            }

          }

        }

      } else if(num_responses[item] > 2) {

        # Intercept updates.
        anl_deriv <- d_categorical("c0",
                                   p[[item]],
                                   etable_item,
                                   theta,
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
                                     theta,
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
                                     theta,
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

        if(!any(item == anchor)){

          p2 <- unlist(p)

          # Intercept DIF updates.
          for(cov in 1:num_predictors) {

            # End routine if only one anchor item is left on each covariate
            # for each item parameter.
            if(is.null(anchor) &
               sum(p2[grep(paste0("c1(.*?)cov",cov),names(p2))] != 0) >
               (num_items - 1) &
               alpha == 1 &&
               (length(final_control$start.values) == 0 || pen > 1) &&
               num_tau >= 10){
              under_identified <- TRUE
              break
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
            if(max_tau) id_max_z <- c(id_max_z,z)
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
               (num_items - 1) &
               alpha == 1 &&
               (length(final_control$start.values) == 0 || pen > 1) &&
               num_tau >= 10){
              under_identified <- TRUE
              break
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
              if(max_tau) id_max_z <- c(id_max_z,z)
              p[[item]][[length(p[[item]])-ncol(pred_data)+cov]] <-
                ifelse(pen_type == "lasso",
                       soft_thresh_cpp(z,alpha,tau_current),
                       firm_thresh_cpp(z,alpha,tau_current,gamma))
            }
          }
        }


        # Gaussian responses.
      } else if(num_responses[item] == 1) {

        # Intercept updates.
        anl_deriv <- d_mu_gaussian("c0",
                                   p[[item]],
                                   etable,
                                   theta,
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
                                     etable_item,
                                     theta,
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
                                      etable_item,
                                      theta,
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
                                          etable_item,
                                          theta,
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
               (num_items - 1) &
               alpha == 1 &&
               (length(final_control$start.values) == 0 || pen > 1) &&
               num_tau >= 10){
              under_identified <- TRUE
              break
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
            z <- p[[item]][c1_parms] - anl_deriv[[1]]/anl_deriv[[2]]
            if(max_tau) id_max_z <- c(id_max_z,z)
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
               (num_items - 1) &
               alpha == 1 &&
               (length(final_control$start.values) == 0 || pen > 1) &&
               num_tau >= 10){
              under_identified <- TRUE
              break
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
              z <- p[[item]][a1_parms] - anl_deriv[[1]]/anl_deriv[[2]]
              if(max_tau) id_max_z <- c(id_max_z,z)
              p[[item]][a1_parms] <- ifelse(pen_type == "lasso",
                                            soft_thresh_cpp(z,alpha,tau_current),
                                            firm_thresh_cpp(z,alpha,tau_current,gamma))
            }
          }
        }

      }

    }

      p_cd_all <- p_cd

      # Update and check for convergence: Calculate the difference
      # in parameter estimates from current to previous.
      eps_cd_all = sqrt(sum((unlist(p_cd_all)-unlist(lastp_cd_all))^2))

      # Update parameter list.
      lastp_cd_all <- p_cd_all

      print(paste0("CD iter: ",iter_cd_all))
      # Update the iteration number.
      iter_cd_all = iter_cd_all + 1

    }

    if(max_tau) {
      id_max_z <- max(abs(id_max_z))
      return(id_max_z)
    } else {
      return(list(p=p_cd_all,
                  under_identified=under_identified))
    }

  }
