#' Maximization step.
#'
#' @param p List of parameters.
#' @param item_data Matrix or data frame of item responses.
#' @param pred_data Matrix or data frame of DIF and/or impact predictors.
#' @param prox_data Vector of observed proxy scores.
#' @param mean_predictors Possibly different matrix of predictors for the mean
#' impact equation.
#' @param var_predictors Possibly different matrix of predictors for the
#' variance impact equation.
#' @param eout E-step output, including matrix for item and impact equations,
#' in addition to theta values (possibly adaptive).
#' @param item_type Optional character value or vector indicating the type of
#' item to be modeled.
#' @param pen_type Character value indicating the penalty function to use.
#' @param tau_current A single numeric value of tau that exists within
#' \code{tau_vec}.
#' @param pen Current penalty index.
#' @param pen.deriv Logical value indicating whether to use the second
#' derivative of the penalized parameter during regularization. The default is
#' TRUE.
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
#' @param max_tau Logical indicating whether to output the minimum tau value
#' needed to remove all DIF from the model.
#'
#' @return a \code{"list"} of estimates obtained from the maximization step using univariate
#' Newton-Raphson
#'
#' @keywords internal
#'
Mstep_simple <-
  function(p,
           item_data,
           pred_data,
           prox_data,
           mean_predictors,
           var_predictors,
           eout,
           item_type,
           pen_type,
           tau_current,
           pen,
           pen.deriv,
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
           max_tau,
           method) {

    # Set under-identified model to FALSE until proven TRUE.
    under_identified <- FALSE

    if(is.null(prox_data)) {
      # Update theta and etable.
      theta <- eout$theta
      etable <- eout$etable
    }


    # Last Mstep.
    if(max_tau) id_max_z <- 0

    # Latent mean impact updates.
    for(cov in 1:ncol(mean_predictors)) {


      if(is.null(prox_data)) {

        anl_deriv <- d_alpha(c(p[[num_items+1]], p[[num_items+2]]),
                             etable,
                             theta,
                             mean_predictors,
                             var_predictors,
                             cov=cov,
                             samp_size,
                             num_items,
                             num_quad)

      } else {

        anl_deriv <- d_alpha_proxy(c(p[[num_items+1]], p[[num_items+2]]),
                                   prox_data,
                                   mean_predictors,
                                   var_predictors,
                                   cov=cov,
                                   samp_size,
                                   num_items)

      }

      p[[num_items+1]][[cov]] <-
        p[[num_items+1]][[cov]] - anl_deriv[[1]]/anl_deriv[[2]]

    }

    # Latent variance impact updates.
    for(cov in 1:ncol(var_predictors)) {

      if(is.null(prox_data)) {

        anl_deriv <- d_phi(c(p[[num_items+1]], p[[num_items+2]]),
                           etable,
                           theta,
                           mean_predictors,
                           var_predictors,
                           cov=cov,
                           samp_size,
                           num_items,
                           num_quad)


      } else {

        anl_deriv <- d_phi_proxy(c(p[[num_items+1]], p[[num_items+2]]),
                                 prox_data,
                                 mean_predictors,
                                 var_predictors,
                                 cov=cov,
                                 samp_size,
                                 num_items)

      }

      p[[num_items+2]][[cov]] <-
        p[[num_items+2]][[cov]] - anl_deriv[[1]]/anl_deriv[[2]]

    }

    for (item in 1:num_items) {


      # Obtain E-tables for each response category.
      if(item_type[item] != "cfa" & is.null(prox_data)) {
        etable_item <- lapply(1:num_responses[item], function(x) etable)
        for(resp in 1:num_responses[item]) {
          etable_item[[resp]][which(
            !(item_data[,item] == resp)), ] <- 0
        }
      }

      if(item_type[item] == "cfa") {

          # Intercept updates.
          if(is.null(prox_data)) {
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
          } else {
            anl_deriv <- d_mu_gaussian_proxy("c0",
                                             p[[item]],
                                             prox_data,
                                             item_data[,item],
                                             pred_data,
                                             cov=NULL,
                                             samp_size)
          }

          p[[item]][[1]] <- p[[item]][[1]] - anl_deriv[[1]]/anl_deriv[[2]]


        # Slope updates.
        if(item_type[item] != "rasch") {

          a0_parms <- grep(paste0("a0_item",item,"_"),names(p[[item]]),fixed=T)


            if(is.null(prox_data)) {
              anl_deriv <- d_mu_gaussian("a0",
                                         p[[item]],
                                         etable,
                                         theta,
                                         item_data[,item],
                                         pred_data,
                                         cov=NULL,
                                         samp_size,
                                         num_items,
                                         num_quad)
            } else {
              anl_deriv <- d_mu_gaussian_proxy("a0",
                                               p[[item]],
                                               prox_data,
                                               item_data[,item],
                                               pred_data,
                                               cov=NULL,
                                               samp_size)
            }

            p[[item]][a0_parms] <- p[[item]][a0_parms] - anl_deriv[[1]]/anl_deriv[[2]]

        } # End Rasch conditional.


        # Residual updates.
        s0_parms <- grep(paste0("s0_item",item,"_"),names(p[[item]]),fixed=T)

          if(is.null(prox_data)) {
            anl_deriv <- d_sigma_gaussian("s0",
                                          p[[item]],
                                          etable,
                                          theta,
                                          item_data[,item],
                                          pred_data,
                                          cov=NULL,
                                          samp_size,
                                          num_items,
                                          num_quad)
          } else {
            anl_deriv <- d_sigma_gaussian_proxy("s0",
                                                p[[item]],
                                                prox_data,
                                                item_data[,item],
                                                pred_data,
                                                cov=NULL,
                                                samp_size)
          }

          p[[item]][s0_parms][[1]] <- p[[item]][s0_parms][[1]] - anl_deriv[[1]]/anl_deriv[[2]]

          if(p[[item]][s0_parms][[1]] < 0) p[[item]][s0_parms][[1]] <- 1


        if(!any(item == anchor)) {


          # Residual DIF updates.
          for(cov in 1:num_predictors) {

            s1_parms <-
              grep(paste0("s1_item",item,"_cov",cov),names(p[[item]]),fixed=T)


              if(is.null(prox_data)) {
                anl_deriv <- d_sigma_gaussian("s1",
                                              p[[item]],
                                              etable,
                                              theta,
                                              item_data[,item],
                                              pred_data,
                                              cov=cov,
                                              samp_size,
                                              num_items,
                                              num_quad)
              } else{
                anl_deriv <- d_sigma_gaussian_proxy("s1",
                                                    p[[item]],
                                                    prox_data,
                                                    item_data[,item],
                                                    pred_data,
                                                    cov=cov,
                                                    samp_size)
              }

              p[[item]][s1_parms][[1]] <- p[[item]][s1_parms][[1]] -
                anl_deriv[[1]]/anl_deriv[[2]]

              if(method == "UNR") {
                p[[item]][s1_parms][[1]] <- p[[item]][s1_parms][[1]]
                break
              }


          } # End looping across covariates.

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
              grep(paste0("c1_item",item,"_cov",cov),names(p[[item]]),fixed=T)

              if(is.null(prox_data)) {
                anl_deriv <- d_mu_gaussian("c1",
                                           p[[item]],
                                           etable,
                                           theta,
                                           item_data[,item],
                                           pred_data,
                                           cov,
                                           samp_size,
                                           num_items,
                                           num_quad)
              } else{
                anl_deriv <- d_mu_gaussian_proxy("c1",
                                                 p[[item]],
                                                 prox_data,
                                                 item_data[,item],
                                                 pred_data,
                                                 cov,
                                                 samp_size)
              }

              z_int <- p[[item]][c1_parms] - anl_deriv[[1]]/anl_deriv[[2]]

              if(max_tau & pen.deriv) {
                id_max_z <- c(id_max_z,z_int*(-anl_deriv[[2]]))
              } else if(max_tau & !pen.deriv) {
                id_max_z <- c(id_max_z,z_int)
              }

              if(!(pen_type == "grp.lasso" || pen_type == "grp.mcp")) {

                if(pen.deriv) {
                  p[[item]][c1_parms][[1]] <-
                    ifelse(pen_type == "lasso",
                           soft_threshold(z_int,alpha,tau_current/-anl_deriv[[2]]),
                           firm_threshold(z_int,alpha,tau_current/-anl_deriv[[2]],gamma))
                } else {
                  p[[item]][c1_parms][[1]] <-
                    ifelse(pen_type == "lasso",
                           soft_threshold(z_int,alpha,tau_current),
                           firm_threshold(z_int,alpha,tau_current,gamma))
                }

              }



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



            if(item_type[item] != "rasch"){

              a1_parms <- grep(paste0("a1_item",item,"_cov",cov),names(p[[item]]),fixed=T)

                if(is.null(prox_data)) {
                  anl_deriv <- d_mu_gaussian("a1",
                                             p[[item]],
                                             etable,
                                             theta,
                                             item_data[,item],
                                             pred_data,
                                             cov,
                                             samp_size,
                                             num_items,
                                             num_quad)
                } else {
                  anl_deriv <- d_mu_gaussian_proxy("a1",
                                                   p[[item]],
                                                   prox_data,
                                                   item_data[,item],
                                                   pred_data,
                                                   cov,
                                                   samp_size)
                }


                z_slp <- p[[item]][a1_parms] - anl_deriv[[1]]/anl_deriv[[2]]

                if(max_tau & pen.deriv) {
                  id_max_z <- c(id_max_z,z_slp*(-anl_deriv[[2]]))
                } else if(max_tau & !pen.deriv) {
                  id_max_z <- c(id_max_z,z_slp)
                }


                if(!(pen_type == "grp.lasso" || pen_type == "grp.mcp")) {

                  if(pen.deriv) {
                    p[[item]][a1_parms][[1]] <-
                      ifelse(pen_type == "lasso",
                             soft_threshold(z_slp,alpha,tau_current/-anl_deriv[[2]]),
                             firm_threshold(z_slp,alpha,tau_current/-anl_deriv[[2]],gamma))
                  } else {
                    p[[item]][a1_parms][[1]] <-
                      ifelse(pen_type == "lasso",
                             soft_threshold(z_slp,alpha,tau_current),
                             firm_threshold(z_slp,alpha,tau_current,gamma))
                  }

                } else if(pen_type == "grp.lasso" || pen_type == "grp.mcp"){


                  grp.update <-
                    if(pen_type == "grp.lasso") {

                      grp_soft_threshold(c(z_int,z_slp),
                                         tau_current)

                    } else if(pen_type == "grp.mcp") {

                      grp_firm_threshold(c(z_int,z_slp),
                                         tau_current,
                                         gamma)

                    }

                    p[[item]][c1_parms][[1]] <- grp.update[[1]]
                    p[[item]][a1_parms][[1]] <- grp.update[[2]]


                }


            } # End Rasch conditional.

          } # End looping across covariates.

        } # End anchor item conditional.



        # Bernoulli responses.
      } else if(item_type[item] == "2pl") {


          # Intercept updates.
          if(is.null(prox_data)) {
            anl_deriv <- d_bernoulli("c0",
                                     p[[item]],
                                     etable_item,
                                     theta,
                                     pred_data,
                                     cov=0,
                                     samp_size,
                                     num_items,
                                     num_quad)
          } else {
            anl_deriv <- d_bernoulli_proxy("c0",
                                           p[[item]],
                                           prox_data,
                                           pred_data,
                                           item_data[,item],
                                           cov=0,
                                           samp_size,
                                           num_items)
          }

          p[[item]][[1]] <- p[[item]][[1]] - anl_deriv[[1]]/anl_deriv[[2]]


        # Slope updates.
        if(item_type[item] != "Rasch") {

            if(is.null(prox_data)) {
              anl_deriv <- d_bernoulli("a0",
                                       p[[item]],
                                       etable_item,
                                       theta,
                                       pred_data,
                                       cov=0,
                                       samp_size,
                                       num_items,
                                       num_quad)
            } else {
              anl_deriv <- d_bernoulli_proxy("a0",
                                             p[[item]],
                                             prox_data,
                                             pred_data,
                                             item_data[,item],
                                             cov=0,
                                             samp_size,
                                             num_items)
            }

            p[[item]][[2]] <- p[[item]][[2]] - anl_deriv[[1]]/anl_deriv[[2]]

        } # End Rasch conditional.


        if(!any(item == anchor)) {

          p2 <- unlist(p)

          # Intercept DIF updates.
          for(cov in 1:num_predictors) {

            # if(pen_type == "grp.lasso" || pen_type == "grp.mcp") {
            #
            #   # End routine if only one anchor item is left on each covariate
            #   # for each item parameter.
            #   if(is.null(anchor) &&
            #      sum(p2[c(grep(paste0("c1(.*?)cov",cov),names(p2)),
            #               grep(paste0("a1(.*?)cov",cov),names(p2)))] != 0) >
            #      (num_items*2 - 1) &&
            #      (length(final_control$start.values) == 0 || pen > 1) &&
            #      num_tau >= 10){
            #     under_identified <- TRUE
            #     break
            #   }
            #
            # }

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


              if(is.null(prox_data)) {
                anl_deriv <- d_bernoulli("c1",
                                         p[[item]],
                                         etable_item,
                                         theta,
                                         pred_data,
                                         cov,
                                         samp_size,
                                         num_items,
                                         num_quad)
              } else {
                anl_deriv <- d_bernoulli_proxy("c1",
                                               p[[item]],
                                               prox_data,
                                               pred_data,
                                               item_data[,item],
                                               cov,
                                               samp_size,
                                               num_items)
              }

              z_int <- p[[item]][[2+cov]] - anl_deriv[[1]]/anl_deriv[[2]]

              if(max_tau & pen.deriv) {
                id_max_z <- c(id_max_z,z_int*(-anl_deriv[[2]]))
              } else if(max_tau & !pen.deriv) {
                id_max_z <- c(id_max_z,z_int)
              }

              if(!(pen_type == "grp.lasso" || pen_type == "grp.mcp")) {

                if(pen.deriv) {
                  p[[item]][[2+cov]] <-
                    ifelse(pen_type == "lasso",
                           soft_threshold(z_int,alpha,tau_current/-anl_deriv[[2]]),
                           firm_threshold(z_int,alpha,tau_current/-anl_deriv[[2]],gamma))
                } else {
                  p[[item]][[2+cov]] <-
                    ifelse(pen_type == "lasso",
                           soft_threshold(z_int,alpha,tau_current),
                           firm_threshold(z_int,alpha,tau_current,gamma))
                }

              }



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


                if(is.null(prox_data)) {
                  anl_deriv <- d_bernoulli("a1",
                                           p[[item]],
                                           etable_item,
                                           theta,
                                           pred_data,
                                           cov,
                                           samp_size,
                                           num_items,
                                           num_quad)
                } else {
                  anl_deriv <- d_bernoulli_proxy("a1",
                                                 p[[item]],
                                                 prox_data,
                                                 pred_data,
                                                 item_data[,item],
                                                 cov,
                                                 samp_size,
                                                 num_items)
                }

                z_slp <- p[[item]][[2+num_predictors+cov]] - anl_deriv[[1]]/anl_deriv[[2]]

                if(max_tau & pen.deriv) {
                  id_max_z <- c(id_max_z,z_slp*(-anl_deriv[[2]]))
                } else if(max_tau & !pen.deriv) {
                  id_max_z <- c(id_max_z,z_slp)
                }

                if(!(pen_type == "grp.lasso" || pen_type == "grp.mcp")) {

                  if(pen.deriv) {
                    p[[item]][[2+num_predictors+cov]] <-
                      ifelse(pen_type == "lasso",
                             soft_threshold(z_slp,alpha,tau_current/-anl_deriv[[2]]),
                             firm_threshold(z_slp,alpha,tau_current/-anl_deriv[[2]],gamma))
                  } else {
                    p[[item]][[2+num_predictors+cov]] <-
                      ifelse(pen_type == "lasso",
                             soft_threshold(z_slp,alpha,tau_current),
                             firm_threshold(z_slp,alpha,tau_current,gamma))
                  }

                } else if(pen_type == "grp.lasso" || pen_type == "grp.mcp") {

                  grp.update <-
                    if(pen_type == "grp.lasso") {

                      grp_soft_threshold(c(z_int,z_slp),
                                         tau_current)

                    } else if(pen_type == "grp.mcp") {

                      grp_firm_threshold(c(z_int,z_slp),
                                         tau_current,
                                         gamma)

                    }

                  p[[item]][[2+cov]] <- grp.update[[1]]
                  p[[item]][[2+num_predictors+cov]] <- grp.update[[2]]

                }


            } # End Rasch conditional.

          } # End looping across covariates.

        } # End anchor item conditional.

        # Categorical.
      } else {


          # Intercept updates.
          if(is.null(prox_data)) {
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
          } else {
            anl_deriv <- d_categorical_proxy("c0",
                                             p[[item]],
                                             prox_data,
                                             pred_data,
                                             item_data[,item],
                                             thr=-1,
                                             cov=-1,
                                             samp_size,
                                             num_responses[[item]],
                                             num_items)
          }
          p[[item]][[1]] <- p[[item]][[1]] - anl_deriv[[1]]/anl_deriv[[2]]

          if(method == "UNR") {
            p[[item]][[1]] <- p[[item]][[1]]
            break
          }


        # Threshold updates.


        for(thr in 2:(num_responses[item]-1)) {


            if(is.null(prox_data)) {
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
            } else {
              anl_deriv <- d_categorical_proxy("c0",
                                               p[[item]],
                                               prox_data,
                                               pred_data,
                                               item_data[,item],
                                               thr=thr,
                                               cov=-1,
                                               samp_size,
                                               num_responses[[item]],
                                               num_items)
            }
            p[[item]][[thr]] <- p[[item]][[thr]] - anl_deriv[[1]]/anl_deriv[[2]]

            if(is.na(p[[item]][[thr]])) {
              p[[item]][[thr]] <-
                mean(sapply(1:num_items, function(x) p[[x]][[thr]]), na.rm = T)

            }



        } # End looping across thresholds.

        # Slope updates.
        if(item_type[item] != "rasch") {


            if(is.null(prox_data)) {
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
            } else {
              anl_deriv <- d_categorical_proxy("a0",
                                               p[[item]],
                                               prox_data,
                                               pred_data,
                                               item_data[,item],
                                               thr=-1,
                                               cov=-1,
                                               samp_size,
                                               num_responses[[item]],
                                               num_items)
            }
            p[[item]][[num_responses[[item]]]] <-
              p[[item]][[num_responses[[item]]]] - anl_deriv[[1]]/anl_deriv[[2]]


        } # End Rasch conditional.

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

              if(is.null(prox_data)) {
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
              } else {
                anl_deriv <- d_categorical_proxy("c1",
                                                 p[[item]],
                                                 prox_data,
                                                 pred_data,
                                                 item_data[,item],
                                                 thr=-1,
                                                 cov,
                                                 samp_size,
                                                 num_responses[[item]],
                                                 num_items)
              }

              z <- p[[item]][[num_responses[[item]]+cov]] - anl_deriv[[1]]/anl_deriv[[2]]

              if(max_tau & pen.deriv) {
                id_max_z <- c(id_max_z,z*(-anl_deriv[[2]]))
              } else if(max_tau & !pen.deriv) {
                id_max_z <- c(id_max_z,z)
              }

              if(pen.deriv) {
                p[[item]][[num_responses[[item]]+cov]] <-
                  ifelse(pen_type == "lasso",
                         soft_threshold(z,alpha,tau_current/-anl_deriv[[2]]),
                         firm_threshold(z,alpha,tau_current/-anl_deriv[[2]],gamma))
              } else {
                p[[item]][[num_responses[[item]]+cov]] <-
                  ifelse(pen_type == "lasso",
                         soft_threshold(z,alpha,tau_current),
                         firm_threshold(z,alpha,tau_current,gamma))
              }


          } # End looping across covariates.

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

            if(item_type[item] != "rasch") {


                if(is.null(prox_data)) {
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
                } else {
                  anl_deriv <- d_categorical_proxy("a1",
                                                   p[[item]],
                                                   prox_data,
                                                   pred_data,
                                                   item_data[,item],
                                                   thr=-1,
                                                   cov,
                                                   samp_size,
                                                   num_responses[[item]],
                                                   num_items)
                }

                z <- p[[item]][[length(p[[item]])-ncol(pred_data)+cov]] -
                  anl_deriv[[1]]/anl_deriv[[2]]

                if(max_tau & pen.deriv) {
                  id_max_z <- c(id_max_z,z*(-anl_deriv[[2]]))
                } else if(max_tau & !pen.deriv) {
                  id_max_z <- c(id_max_z,z)
                }


                if(pen.deriv) {
                  p[[item]][[length(p[[item]])-ncol(pred_data)+cov]] <-
                    ifelse(pen_type == "lasso",
                           soft_threshold(z,alpha,tau_current/-anl_deriv[[2]]),
                           firm_threshold(z,alpha,tau_current/-anl_deriv[[2]],gamma))
                } else {
                  p[[item]][[length(p[[item]])-ncol(pred_data)+cov]] <-
                    ifelse(pen_type == "lasso",
                           soft_threshold(z,alpha,tau_current),
                           firm_threshold(z,alpha,tau_current,gamma))
                }



            } # End Rasch conditional.

          } # End looping across covariates.

        } # End anchor item condtional.


      } # End item type conditional.

    } # End looping through items.




    if(max_tau) {

      id_max_z <- max(abs(unlist(id_max_z)))

      return(id_max_z)

    } else {

        return(list(p=p,
                    under_identified=under_identified))
    }

  }
