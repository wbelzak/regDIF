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
#' @param method Character value indicating the type of optimization method. Options include "MNR",
#' "UNR", and "CD"
#'
#' @return a \code{"list"} of estimates obtained from the maximization step using multivariate
#' Newton-Raphson
#'
#' @importFrom foreach %dopar%
#'
#' @keywords internal
#'
Mstep <-
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

    # CD Maximization and print settings.
    lastp_cd_all <- p_cd <- p
    eps_cd_all <- Inf
    iter_cd_all <- 1

    # Loop until convergence or maximum number of iterations.
    while(eps_cd_all > final_control$tol){


    # Latent impact updates.
    if(method == "MNR") {

      if(is.null(prox_data)) {
        anl_deriv_impact <- d_impact_block(p[[num_items+1]],
                                           p[[num_items+2]],
                                           etable,
                                           theta,
                                           mean_predictors,
                                           var_predictors,
                                           samp_size,
                                           num_items,
                                           num_quad,
                                           num_predictors)
      } else {
        anl_deriv_impact <- d_impact_block_proxy(p[[num_items+1]],
                                                 p[[num_items+2]],
                                                 prox_data,
                                                 mean_predictors,
                                                 var_predictors,
                                                 samp_size,
                                                 num_items,
                                                 num_quad,
                                                 num_predictors)
      }


      inv_hess_impact <- solve(anl_deriv_impact[[2]])
      inv_hess_impact_diag <- -diag(inv_hess_impact)
      m <- c(p[[num_items+1]],p[[num_items+2]]) -
        inv_hess_impact %*% anl_deriv_impact[[1]]
      names(m) <- names(c(p[[num_items+1]],p[[num_items+2]]))


      p[[num_items+1]] <- m[1:ncol(mean_predictors)]
      p[[num_items+2]] <- m[(ncol(mean_predictors)+1):length(m)]


      inv_hess_diag <- vector('list',num_items+2)

    } else { # method == "UNR" || method == "CD"

      # Latent mean impact updates.
      for(cov in 1:ncol(mean_predictors)) {

        lastp_cd <- p_cd
        eps_cd <- Inf
        iter_cd <- 1

        # Loop until convergence or maximum number of iterations.
        while(eps_cd > final_control$tol){

          if(is.null(prox_data)) {

            anl_deriv <- d_alpha(c(p_cd[[num_items+1]],p_cd[[num_items+2]]),
                                 etable,
                                 theta,
                                 mean_predictors,
                                 var_predictors,
                                 cov=cov,
                                 samp_size,
                                 num_items,
                                 num_quad)

          } else {

            anl_deriv <- d_alpha_proxy(c(p_cd[[num_items+1]],p_cd[[num_items+2]]),
                                       prox_data,
                                       mean_predictors,
                                       var_predictors,
                                       cov=cov,
                                       samp_size,
                                       num_items)

          }

          p_cd[[num_items+1]][[cov]] <-
            p_cd[[num_items+1]][[cov]] - anl_deriv[[1]]/anl_deriv[[2]]

          if(method == "UNR") {
            p[[num_items+1]][[cov]] <- p_cd[[num_items+1]][[cov]]
            break
          }

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

          if(is.null(prox_data)) {

            anl_deriv <- d_phi(c(p_cd[[num_items+1]],p_cd[[num_items+2]]),
                               etable,
                               theta,
                               mean_predictors,
                               var_predictors,
                               cov=cov,
                               samp_size,
                               num_items,
                               num_quad)


          } else {

            anl_deriv <- d_phi_proxy(c(p_cd[[num_items+1]],p_cd[[num_items+2]]),
                                     prox_data,
                                     mean_predictors,
                                     var_predictors,
                                     cov=cov,
                                     samp_size,
                                     num_items)

          }

          p_cd[[num_items+2]][[cov]] <-
            p_cd[[num_items+2]][[cov]] - anl_deriv[[1]]/anl_deriv[[2]]


          if(method == "UNR") {
            p[[num_items+2]][[cov]] <- p_cd[[num_items+2]][[cov]]
            break
          }

          # Update and check for convergence: Calculate the difference
          # in parameter estimates from current to previous.
          eps_cd = sqrt(sum((unlist(p_cd)-unlist(lastp_cd))^2))

          # Update parameter list.
          lastp_cd <- p_cd

          # Update the iteration number.
          iter_cd = iter_cd + 1
        }

      }

    } # End optimization method.


    # Item response updates.
    if(method == "MNR") {

        for (item in 1:num_items) {

          # Gaussian responses
          if (item_type[item] == "cfa") {


            if(is.null(prox_data)) {
              anl_deriv_item <- d_gaussian_itemblock(p[[item]],
                                                     etable,
                                                     theta,
                                                     item_data[,item],
                                                     pred_data,
                                                     samp_size,
                                                     num_items,
                                                     num_quad,
                                                     num_predictors)
            } else {
              anl_deriv_item <- d_gaussian_itemblock_proxy(p[[item]],
                                                           prox_data,
                                                           item_data[,item],
                                                           pred_data,
                                                           samp_size,
                                                           num_items,
                                                           num_quad,
                                                           num_predictors)
            }


            inv_hess_item <- solve(anl_deriv_item[[2]])
            inv_hess_diag[[item]] <- -diag(inv_hess_item)
            z <- p[[item]] - inv_hess_item %*% anl_deriv_item[[1]]

            # c0 update.
            p[[item]][[1]] <- z[1]

            # a0 update.
            if(item_type[item] != "rasch") p[[item]][[2]] <- z[2]

            # s0 update.
            p[[item]][[3+num_predictors*2]] <- z[3+num_predictors*2]
            if(p[[item]][[3+num_predictors*2]] < 0) {
              p[[item]][[3+num_predictors*2]] <- z[3+num_predictors*2] <- 1
            }

            # Don't update DIF estimates if anchor item.
            if(any(item == anchor)) next

            p2 <- unlist(p)

            for(cov in 1:num_predictors) {

              # If penalty type is a group function.
              if(pen_type == "grp.lasso" ||
                 pen_type == "grp.mcp") {

                # End routine if only one anchor item is left on each covariate
                # for each item parameter.
                if(is.null(anchor) &&
                   sum(p2[c(grep(paste0("c1(.*?)cov",cov),names(p2)),
                            grep(paste0("a1(.*?)cov",cov),names(p2)))] != 0) >
                   (num_items*2 - 1) &&
                   (length(final_control$start.values) == 0 || pen > 1) &&
                   num_tau >= 10){
                  under_identified <- TRUE
                  break
                }

                if(max_tau) {
                  id_max_z <- c(id_max_z,
                                z[2+cov],
                                z[2+num_predictors+cov])
                }


                # group update.
                grp.update <-
                  if(pen_type == "grp.lasso") {
                    grp_soft_threshold(z[c(2+cov,
                                           2+num_predictors+cov)],
                                       tau_current)
                  } else if(pen_type == "grp.mcp") {
                    grp_firm_threshold(z[c(2+cov,
                                           2+num_predictors+cov)],
                                       tau_current,
                                       gamma)
                  }

                # c1 updates.
                p[[item]][[2+cov]] <- grp.update[[1]]
                # a1 updates.
                p[[item]][[2+num_predictors+cov]] <- grp.update[[2]]

                # s1 updates.
                p[[item]][[2+num_predictors*2+cov]] <- z[2+num_predictors*2+cov]

                next

              } # End group penalty conditional.

              # End routine if only one anchor item is left on each covariate
              # for each item parameter.
              if(is.null(anchor) &&
                 sum(p2[grep(paste0("c1(.*?)cov",cov),names(p2))] != 0) >
                 (num_items - 1) &&
                 alpha == 1 &&
                 (length(final_control$start.values) == 0 || pen > 1) &&
                 num_tau >= 10){
                under_identified <- TRUE
                break
              }

              if(max_tau) {
                id_max_z <- c(id_max_z,
                              z[2+cov],
                              z[2+num_predictors+cov])
              }

              # c1 updates.
              p[[item]][[2+cov]] <-
                if(pen_type == "lasso") {
                  soft_threshold(z[2+cov],
                                 alpha,
                                 tau_current)
                } else if(pen_type == "mcp") {
                  firm_threshold(z[2+cov],
                                 alpha,
                                 tau_current,
                                 gamma)
                }

              # s1 updates.
              p[[item]][[2+num_predictors*2+cov]] <- z[2+num_predictors*2+cov]


              # a1 updates.
              if(item_type[item] == "rasch") next

              # End routine if only one anchor item is left on each covariate
              # for each item parameter.
              if(is.null(anchor) &&
                 sum(p2[grep(paste0("a1(.*?)cov",cov),names(p2))] != 0) >
                 (num_items - 1) &&
                 alpha == 1 &&
                 (length(final_control$start.values) == 0 || pen > 1) &&
                 num_tau >= 10){
                under_identified <- TRUE
                break
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

            } # End looping across covariates.

            if(under_identified) break


          } else if(item_type[item] == "2pl") {

            if(is.null(prox_data)) {
              anl_deriv_item <- d_bernoulli_itemblock(p[[item]],
                                                      etable,
                                                      theta,
                                                      pred_data,
                                                      item_data[,item],
                                                      samp_size,
                                                      num_items,
                                                      num_predictors,
                                                      num_quad)
            } else {
              anl_deriv_item <- d_bernoulli_itemblock_proxy(p[[item]],
                                                            pred_data,
                                                            item_data[,item],
                                                            prox_data,
                                                            samp_size,
                                                            num_items,
                                                            num_predictors,
                                                            num_quad)
            }


            inv_hess_item <- solve(anl_deriv_item[[2]])
            inv_hess_diag[[item]] <- -diag(inv_hess_item)
            z <- p[[item]] - inv_hess_item %*% anl_deriv_item[[1]]

            # c0 update.
            p[[item]][[1]] <- z[1]

            # a0 update.
            if(item_type[item] != "rasch") p[[item]][[2]] <- z[2]

            # Don't update DIF estimates if anchor item.
            if(any(item == anchor)) next

            p2 <- unlist(p)

            for(cov in 1:num_predictors) {

              # If penalty type is a group function.
              if(pen_type == "grp.lasso" ||
                 pen_type == "grp.mcp") {

                # End routine if only one anchor item is left on each covariate
                # for each item parameter.
                if(is.null(anchor) &&
                   sum(p2[c(grep(paste0("c1(.*?)cov",cov),names(p2)),
                            grep(paste0("a1(.*?)cov",cov),names(p2)))] != 0) >
                   (num_items*2 - 1) &&
                   (length(final_control$start.values) == 0 || pen > 1) &&
                   num_tau >= 10){
                  under_identified <- TRUE
                  break
                }

                if(max_tau) {
                  id_max_z <- c(id_max_z,
                                z[2+cov],
                                z[2+num_predictors+cov])
                }


                # group update.
                grp.update <-
                  if(pen_type == "grp.lasso") {
                    grp_soft_threshold(z[c(2+cov,
                                           2+num_predictors+cov)],
                                       tau_current)
                  } else if(pen_type == "grp.mcp") {
                    grp_firm_threshold(z[c(2+cov,
                                           2+num_predictors+cov)],
                                       tau_current,
                                       gamma)
                  }

                # c1 updates.
                p[[item]][[2+cov]] <- grp.update[[1]]
                # a1 updates.
                p[[item]][[2+num_predictors+cov]] <- grp.update[[2]]

                next

              } # End group penalty conditional.

              # End routine if only one anchor item is left on each covariate
              # for each item parameter.
              if(is.null(anchor) &&
                 sum(p2[grep(paste0("c1(.*?)cov",cov),names(p2))] != 0) >
                 (num_items - 1) &&
                 alpha == 1 &&
                 (length(final_control$start.values) == 0 || pen > 1) &&
                 num_tau >= 10){
                under_identified <- TRUE
                break
              }

              if(max_tau) {
                id_max_z <- c(id_max_z,
                              z[2+cov],
                              z[2+num_predictors+cov])
              }

              # c1 updates.
              p[[item]][[2+cov]] <-
                if(pen_type == "lasso") {
                  soft_threshold(z[2+cov],
                                 alpha,
                                 tau_current)
                } else if(pen_type == "mcp") {
                  firm_threshold(z[2+cov],
                                 alpha,
                                 tau_current,
                                 gamma)
                }


              # a1 updates.
              if(item_type[item] == "rasch") next

              # End routine if only one anchor item is left on each covariate
              # for each item parameter.
              if(is.null(anchor) &&
                 sum(p2[grep(paste0("a1(.*?)cov",cov),names(p2))] != 0) >
                 (num_items - 1) &&
                 alpha == 1 &&
                 (length(final_control$start.values) == 0 || pen > 1) &&
                 num_tau >= 10){
                under_identified <- TRUE
                break
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

            } # End looping across covariates.

            if(under_identified) break

          } # End item type conditional.

        } # End looping across items.


    } else {  # method == "UNR" || method == "CD"

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

          # CD Maximization and print settings.
          lastp_cd <- p_cd
          eps_cd <- Inf
          iter_cd <- 1

          # Loop until convergence or maximum number of iterations.
          while(eps_cd > final_control$tol){

            # Intercept updates.
            if(is.null(prox_data)) {
              anl_deriv <- d_mu_gaussian("c0",
                                         p_cd[[item]],
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
                                               p_cd[[item]],
                                               prox_data,
                                               item_data[,item],
                                               pred_data,
                                               cov=NULL,
                                               samp_size)
            }

            p_cd[[item]][[1]] <- p_cd[[item]][[1]] - anl_deriv[[1]]/anl_deriv[[2]]

            if(method == "UNR") {
              p[[item]][[1]] <- p_cd[[item]][[1]]
              break
            }

            # Update and check for convergence: Calculate the difference
            # in parameter estimates from current to previous.
            eps_cd = sqrt(sum((unlist(p_cd)-unlist(lastp_cd))^2))

            # Update parameter list.
            lastp_cd <- p_cd

            # Update the iteration number.
            iter_cd = iter_cd + 1

          } # End coordinate descent looping.

          # Slope updates.
          if(item_type[item] != "rasch") {

            a0_parms <- grep(paste0("a0_item",item,"_"),names(p_cd[[item]]),fixed=T)

            # CD Maximization and print settings.
            lastp_cd <- p_cd
            eps_cd <- Inf
            iter_cd <- 1

            # Loop until convergence or maximum number of iterations.
            while(eps_cd > final_control$tol){

              if(is.null(prox_data)) {
                anl_deriv <- d_mu_gaussian("a0",
                                           p_cd[[item]],
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
                                                 p_cd[[item]],
                                                 prox_data,
                                                 item_data[,item],
                                                 pred_data,
                                                 cov=NULL,
                                                 samp_size)
              }

              p_cd[[item]][a0_parms] <- p_cd[[item]][a0_parms] - anl_deriv[[1]]/anl_deriv[[2]]

              if(method == "UNR") {
                p[[item]][a0_parms] <- p_cd[[item]][a0_parms]
                break
              }

              # Update and check for convergence: Calculate the difference
              # in parameter estimates from current to previous.
              eps_cd = sqrt(sum((unlist(p_cd)-unlist(lastp_cd))^2))

              # Update parameter list.
              lastp_cd <- p_cd

              # Update the iteration number.
              iter_cd = iter_cd + 1

            } # End coordinate descent looping.

          } # End Rasch conditional.


          # Residual updates.
          s0_parms <- grep(paste0("s0_item",item,"_"),names(p_cd[[item]]),fixed=T)

          # CD Maximization and print settings.
          lastp_cd <- p_cd
          eps_cd <- Inf
          iter_cd <- 1

          # Loop until convergence or maximum number of iterations.
          while(eps_cd > final_control$tol){

            if(is.null(prox_data)) {
              anl_deriv <- d_sigma_gaussian("s0",
                                            p_cd[[item]],
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
                                                  p_cd[[item]],
                                                  prox_data,
                                                  item_data[,item],
                                                  pred_data,
                                                  cov=NULL,
                                                  samp_size)
            }

            p_cd[[item]][s0_parms][[1]] <- p_cd[[item]][s0_parms][[1]] -
              anl_deriv[[1]]/anl_deriv[[2]]
            if(p_cd[[item]][s0_parms][[1]] < 0) p_cd[[item]][s0_parms][[1]] <- 1

            if(method == "UNR") {
              p[[item]][s0_parms][[1]] <- p_cd[[item]][s0_parms][[1]]
              break
            }

            # Update and check for convergence: Calculate the difference
            # in parameter estimates from current to previous.
            eps_cd = sqrt(sum((unlist(p_cd)-unlist(lastp_cd))^2))

            # Update parameter list.
            lastp_cd <- p_cd

            # Update the iteration number.
            iter_cd = iter_cd + 1

          } # End coordinate descent looping.

          if(!any(item == anchor)) {


            # Residual DIF updates.
            for(cov in 1:num_predictors) {

              s1_parms <-
                grep(paste0("s1_item",item,"_cov",cov),names(p_cd[[item]]),fixed=T)

              # CD Maximization and print settings.
              lastp_cd <- p_cd
              eps_cd <- Inf
              iter_cd <- 1

              # Loop until convergence or maximum number of iterations.
              while(eps_cd > final_control$tol){

                if(is.null(prox_data)) {
                  anl_deriv <- d_sigma_gaussian("s1",
                                                p_cd[[item]],
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
                                                      p_cd[[item]],
                                                      prox_data,
                                                      item_data[,item],
                                                      pred_data,
                                                      cov=cov,
                                                      samp_size)
                }

                p_cd[[item]][s1_parms][[1]] <- p_cd[[item]][s1_parms][[1]] -
                  anl_deriv[[1]]/anl_deriv[[2]]

                if(method == "UNR") {
                  p[[item]][s1_parms][[1]] <- p_cd[[item]][s1_parms][[1]]
                  break
                }

                # Update and check for convergence: Calculate the difference
                # in parameter estimates from current to previous.
                eps_cd = sqrt(sum((unlist(p_cd)-unlist(lastp_cd))^2))

                # Update parameter list.
                lastp_cd <- p_cd

                # Update the iteration number.
                iter_cd = iter_cd + 1

              } # End coordinate descent looping.

            } # End looping across covariates.

            p2_cd <- unlist(p_cd)

            # Intercept DIF updates.
            for(cov in 1:num_predictors){

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

              c1_parms <-
                grep(paste0("c1_item",item,"_cov",cov),names(p_cd[[item]]),fixed=T)

              # CD Maximization and print settings.
              lastp_cd <- p_cd
              eps_cd <- Inf
              iter_cd <- 1

              # Loop until convergence or maximum number of iterations.
              while(eps_cd > final_control$tol){

                if(is.null(prox_data)) {
                  anl_deriv <- d_mu_gaussian("c1",
                                             p_cd[[item]],
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
                                                   p_cd[[item]],
                                                   prox_data,
                                                   item_data[,item],
                                                   pred_data,
                                                   cov,
                                                   samp_size)
                }

                z <- p_cd[[item]][c1_parms] - anl_deriv[[1]]/anl_deriv[[2]]
                if(max_tau) id_max_z <- c(id_max_z,z)
                p_cd[[item]][c1_parms][[1]] <-
                  ifelse(pen_type == "lasso",
                         soft_threshold(z,alpha,tau_current),
                         firm_threshold(z,alpha,tau_current,gamma))

                if(method == "UNR") {
                  p[[item]][c1_parms] <- p_cd[[item]][c1_parms]
                  break
                }

                # Update and check for convergence: Calculate the difference
                # in parameter estimates from current to previous.
                eps_cd = sqrt(sum((unlist(p_cd)-unlist(lastp_cd))^2))

                # Update parameter list.
                lastp_cd <- p_cd

                # Update the iteration number.
                iter_cd = iter_cd + 1

              } # End coordinate descent looping.

            } # End looping across covariates.

            # Slope DIF updates.
            for(cov in 1:num_predictors){

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

              if(item_type[item] != "rasch"){
                a1_parms <-
                  grep(paste0("a1_item",item,"_cov",cov),names(p_cd[[item]]),fixed=T)

                # CD Maximization and print settings.
                lastp_cd <- p_cd
                eps_cd <- Inf
                iter_cd <- 1

                # Loop until convergence or maximum number of iterations.
                while(eps_cd > final_control$tol){

                  if(is.null(prox_data)) {
                    anl_deriv <- d_mu_gaussian("a1",
                                               p_cd[[item]],
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
                                                     p_cd[[item]],
                                                     prox_data,
                                                     item_data[,item],
                                                     pred_data,
                                                     cov,
                                                     samp_size)
                  }

                  z <- p_cd[[item]][a1_parms] - anl_deriv[[1]]/anl_deriv[[2]]
                  if(max_tau) id_max_z <- c(id_max_z,z)
                  p_cd[[item]][a1_parms][[1]] <- ifelse(pen_type == "lasso",
                                                     soft_threshold(z,alpha,tau_current),
                                                     firm_threshold(z,alpha,tau_current,gamma))



                  if(method == "UNR") {
                    p_cd[[item]][c1_parms] <- p_cd[[item]][a1_parms]
                    break
                  }

                  # Update and check for convergence: Calculate the difference
                  # in parameter estimates from current to previous.
                  eps_cd = sqrt(sum((unlist(p_cd)-unlist(lastp_cd))^2))

                  # Update parameter list.
                  lastp_cd <- p_cd

                  # Update the iteration number.
                  iter_cd = iter_cd + 1

                } # End coordinate descent looping.


              } # End Rasch conditional.

            } # End looping across covariates.

          } # End anchor item conditional.



          # Bernoulli responses.
        } else if(item_type[item] == "2pl") {

          # CD Maximization and print settings.
          lastp_cd <- p_cd
          eps_cd <- Inf
          iter_cd <- 1

          # Loop until convergence or maximum number of iterations.
          while(eps_cd > final_control$tol){

            # Intercept updates.
            if(is.null(prox_data)) {
              anl_deriv <- d_bernoulli("c0",
                                       p_cd[[item]],
                                       etable_item,
                                       theta,
                                       pred_data,
                                       cov=0,
                                       samp_size,
                                       num_items,
                                       num_quad)
            } else {
              anl_deriv <- d_bernoulli_proxy("c0",
                                             p_cd[[item]],
                                             prox_data,
                                             pred_data,
                                             item_data[,item],
                                             cov=0,
                                             samp_size,
                                             num_items)
            }

            p_cd[[item]][[1]] <- p_cd[[item]][[1]] - anl_deriv[[1]]/anl_deriv[[2]]

            if(method == "UNR") {
              p[[item]][[1]] <- p_cd[[item]][[1]]
              break
            }

            # Update and check for convergence: Calculate the difference
            # in parameter estimates from current to previous.
            eps_cd = sqrt(sum((unlist(p_cd)-unlist(lastp_cd))^2))

            # Update parameter list.
            lastp_cd <- p_cd

            # Update the iteration number.
            iter_cd = iter_cd + 1

          } # End coordinate descent looping.

          # Slope updates.
          if(item_type[item] != "Rasch") {

            # CD Maximization and print settings.
            lastp_cd <- p_cd
            eps_cd <- Inf
            iter_cd <- 1

            # Loop until convergence or maximum number of iterations.
            while(eps_cd > final_control$tol){

              if(is.null(prox_data)) {
                anl_deriv <- d_bernoulli("a0",
                                         p_cd[[item]],
                                         etable_item,
                                         theta,
                                         pred_data,
                                         cov=0,
                                         samp_size,
                                         num_items,
                                         num_quad)
              } else {
                anl_deriv <- d_bernoulli_proxy("a0",
                                               p_cd[[item]],
                                               prox_data,
                                               pred_data,
                                               item_data[,item],
                                               cov=0,
                                               samp_size,
                                               num_items)
              }

              p_cd[[item]][[2]] <- p_cd[[item]][[2]] - anl_deriv[[1]]/anl_deriv[[2]]

              if(method == "UNR") {
                p[[item]][[2]] <- p_cd[[item]][[2]]
                break
              }

              # Update and check for convergence: Calculate the difference
              # in parameter estimates from current to previous.
              eps_cd = sqrt(sum((unlist(p_cd)-unlist(lastp_cd))^2))

              # Update parameter list.
              lastp_cd <- p_cd

              # Update the iteration number.
              iter_cd = iter_cd + 1

            } # End coordinate descent looping.

          } # End Rasch conditional.


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

                if(is.null(prox_data)) {
                  anl_deriv <- d_bernoulli("c1",
                                           p_cd[[item]],
                                           etable_item,
                                           theta,
                                           pred_data,
                                           cov,
                                           samp_size,
                                           num_items,
                                           num_quad)
                } else {
                  anl_deriv <- d_bernoulli_proxy("c1",
                                                 p_cd[[item]],
                                                 prox_data,
                                                 pred_data,
                                                 item_data[,item],
                                                 cov,
                                                 samp_size,
                                                 num_items)
                }

                z <- p_cd[[item]][[2+cov]] - anl_deriv[[1]]/anl_deriv[[2]]
                if(max_tau) id_max_z <- c(id_max_z,z)
                p_cd[[item]][[2+cov]] <- ifelse(pen_type == "lasso",
                                                soft_threshold(z,alpha,tau_current),
                                                firm_threshold(z,alpha,tau_current,gamma))

                if(method == "UNR") {
                  p[[item]][[2+cov]] <- p_cd[[item]][[2+cov]]
                  break
                }

                # Update and check for convergence: Calculate the difference
                # in parameter estimates from current to previous.
                eps_cd = sqrt(sum((unlist(p_cd)-unlist(lastp_cd))^2))

                # Update parameter list.
                lastp_cd <- p_cd

                # Update the iteration number.
                iter_cd = iter_cd + 1

              } # End coordinate descent looping.

            } # End looping across covariates.

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

                  if(is.null(prox_data)) {
                    anl_deriv <- d_bernoulli("a1",
                                             p_cd[[item]],
                                             etable_item,
                                             theta,
                                             pred_data,
                                             cov,
                                             samp_size,
                                             num_items,
                                             num_quad)
                  } else {
                    anl_deriv <- d_bernoulli_proxy("a1",
                                                   p_cd[[item]],
                                                   prox_data,
                                                   pred_data,
                                                   item_data[,item],
                                                   cov,
                                                   samp_size,
                                                   num_items)
                  }

                  z <- p_cd[[item]][[2+num_predictors+cov]] -
                    anl_deriv[[1]]/anl_deriv[[2]]
                  if(max_tau) id_max_z <- c(id_max_z,z)
                  p_cd[[item]][[2+num_predictors+cov]] <-
                    ifelse(pen_type == "lasso",
                           soft_threshold(z,alpha,tau_current),
                           firm_threshold(z,alpha,tau_current,gamma))

                  if(method == "UNR") {
                    p[[item]][[2+num_predictors+cov]] <- p_cd[[item]][[2+num_predictors+cov]]
                    break
                  }

                  # Update and check for convergence: Calculate the difference
                  # in parameter estimates from current to previous.
                  eps_cd = sqrt(sum((unlist(p_cd)-unlist(lastp_cd))^2))

                  # Update parameter list.
                  lastp_cd <- p_cd

                  # Update the iteration number.
                  iter_cd = iter_cd + 1

                } # End coordinate descent looping.

              } # End Rasch conditional.

            } # End looping across covariates.

          } # End anchor item conditional.

          # Categorical.
        } else {

          # CD Maximization and print settings.
          lastp_cd <- p_cd
          eps_cd <- Inf
          iter_cd <- 1

          # Loop until convergence or maximum number of iterations.
          while(eps_cd > final_control$tol){

            # Intercept updates.
            if(is.null(prox_data)) {
              anl_deriv <- d_categorical("c0",
                                         p_cd[[item]],
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
                                               p_cd[[item]],
                                               prox_data,
                                               pred_data,
                                               item_data[,item],
                                               thr=-1,
                                               cov=-1,
                                               samp_size,
                                               num_responses[[item]],
                                               num_items)
            }
            p_cd[[item]][[1]] <- p_cd[[item]][[1]] - anl_deriv[[1]]/anl_deriv[[2]]

            if(method == "UNR") {
              p[[item]][[1]] <- p_cd[[item]][[1]]
              break
            }

            # Update and check for convergence: Calculate the difference
            # in parameter estimates from current to previous.
            eps_cd = sqrt(sum((unlist(p_cd)-unlist(lastp_cd))^2))

            # Update parameter list.
            lastp_cd <- p_cd

            # Update the iteration number.
            iter_cd = iter_cd + 1

          } # End coordinate descent looping.


          # Threshold updates.


          for(thr in 2:(num_responses[item]-1)) {
            # avg_thr <- mean(sapply(1:num_items, function(x) p_cd[[x]][[thr]]))

            # CD Maximization and print settings.
            lastp_cd <- p_cd
            eps_cd <- Inf
            iter_cd <- 1

            # Loop until convergence or maximum number of iterations.
            while(eps_cd > final_control$tol){


              if(is.null(prox_data)) {
                anl_deriv <- d_categorical("c0",
                                           p_cd[[item]],
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
                                                 p_cd[[item]],
                                                 prox_data,
                                                 pred_data,
                                                 item_data[,item],
                                                 thr=thr,
                                                 cov=-1,
                                                 samp_size,
                                                 num_responses[[item]],
                                                 num_items)
              }
              p_cd[[item]][[thr]] <- p_cd[[item]][[thr]] - anl_deriv[[1]]/anl_deriv[[2]]

              if(is.na(p_cd[[item]][[thr]])) {
                p_cd[[item]][[thr]] <-
                  mean(sapply(1:num_items, function(x) p_cd[[x]][[thr]]), na.rm = T)
                if(method == "UNR") {
                  p[[item]][[thr]] <- p_cd[[item]][[thr]]
                }
                break
              }


              # if(abs(p_cd[[item]][[thr]]) > 10*avg_thr) {
              #   p_cd[[item]][[thr]] <- avg_thr
              #   warning(paste0("Problem with estimating threshold ", thr, " for item ",
              #                  item,".\n",
              #                  "Setting threshold ", thr, " for item ", item, " to average ",
              #                  "of all other item thresholds."),
              #           call. = FALSE, immediate. = TRUE)
              #   if(method == "UNR") {
              #     p[[item]][[thr]] <- p_cd[[item]][[thr]]
              #   }
              #   break
              # }

              if(method == "UNR") {
                p[[item]][[thr]] <- p_cd[[item]][[thr]]
                break
              }

              # print(p_cd[[item]][3])
              # p_cd[[item]][3]=.25

              # Update and check for convergence: Calculate the difference
              # in parameter estimates from current to previous.
              eps_cd = sqrt(sum((unlist(p_cd)-unlist(lastp_cd))^2))

              # Update parameter list.
              lastp_cd <- p_cd

              # Update the iteration number.
              iter_cd = iter_cd + 1

            } # End coordinate descent looping.


          } # End looping across thresholds.

          # Slope updates.
          if(item_type[item] != "rasch") {

            # CD Maximization and print settings.
            lastp_cd <- p_cd
            eps_cd <- Inf
            iter_cd <- 1

            # Loop until convergence or maximum number of iterations.
            while(eps_cd > final_control$tol){

              if(is.null(prox_data)) {
                anl_deriv <- d_categorical("a0",
                                           p_cd[[item]],
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
                                                 p_cd[[item]],
                                                 prox_data,
                                                 pred_data,
                                                 item_data[,item],
                                                 thr=-1,
                                                 cov=-1,
                                                 samp_size,
                                                 num_responses[[item]],
                                                 num_items)
              }
              p_cd[[item]][[num_responses[[item]]]] <-
                p_cd[[item]][[num_responses[[item]]]] - anl_deriv[[1]]/anl_deriv[[2]]

              if(method == "UNR") {
                p[[item]][[num_responses[[item]]]] <- p_cd[[item]][[num_responses[[item]]]]
                break
              }

              # Update and check for convergence: Calculate the difference
              # in parameter estimates from current to previous.
              eps_cd = sqrt(sum((unlist(p_cd)-unlist(lastp_cd))^2))

              # Update parameter list.
              lastp_cd <- p_cd

              # Update the iteration number.
              iter_cd = iter_cd + 1

            } # End coordinate descent looping.

          } # End Rasch conditional.

          if(!any(item == anchor)){

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

                if(is.null(prox_data)) {
                  anl_deriv <- d_categorical("c1",
                                             p_cd[[item]],
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
                                                   p_cd[[item]],
                                                   prox_data,
                                                   pred_data,
                                                   item_data[,item],
                                                   thr=-1,
                                                   cov,
                                                   samp_size,
                                                   num_responses[[item]],
                                                   num_items)
                }

                z <- p_cd[[item]][[num_responses[[item]]+cov]] -
                  anl_deriv[[1]]/anl_deriv[[2]]
                if(max_tau) id_max_z <- c(id_max_z,z)
                p_cd[[item]][[num_responses[[item]]+cov]] <-
                  ifelse(pen_type == "lasso",
                         soft_threshold(z,alpha,tau_current),
                         firm_threshold(z,alpha,tau_current,gamma))

                if(method == "UNR") {
                  p[[item]][[num_responses[[item]]+cov]] <-
                    p_cd[[item]][[num_responses[[item]]+cov]]
                  break
                }

                # Update and check for convergence: Calculate the difference
                # in parameter estimates from current to previous.
                eps_cd = sqrt(sum((unlist(p_cd)-unlist(lastp_cd))^2))

                # Update parameter list.
                lastp_cd <- p_cd

                # Update the iteration number.
                iter_cd = iter_cd + 1

              } # End coordinate descent looping.



            } # End looping across covariates.

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

              if(item_type[item] != "rasch") {

                # CD Maximization and print settings.
                lastp_cd <- p_cd
                eps_cd <- Inf
                iter_cd <- 1

                # Loop until convergence or maximum number of iterations.
                while(eps_cd > final_control$tol){

                  if(is.null(prox_data)) {
                    anl_deriv <- d_categorical("a1",
                                               p_cd[[item]],
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
                                                     p_cd[[item]],
                                                     prox_data,
                                                     pred_data,
                                                     item_data[,item],
                                                     thr=-1,
                                                     cov,
                                                     samp_size,
                                                     num_responses[[item]],
                                                     num_items)
                  }

                  z <- p_cd[[item]][[length(p[[item]])-ncol(pred_data)+cov]] -
                    anl_deriv[[1]]/anl_deriv[[2]]
                  if(max_tau) id_max_z <- c(id_max_z,z)
                  p_cd[[item]][[length(p[[item]])-ncol(pred_data)+cov]] <-
                    ifelse(pen_type == "lasso",
                           soft_threshold(z,alpha,tau_current),
                           firm_threshold(z,alpha,tau_current,gamma))

                  if(method == "UNR") {
                    p[[item]][[length(p[[item]])-ncol(pred_data)+cov]] <-
                      p_cd[[item]][[length(p[[item]])-ncol(pred_data)+cov]]
                    break
                  }

                  # Update and check for convergence: Calculate the difference
                  # in parameter estimates from current to previous.
                  eps_cd = sqrt(sum((unlist(p_cd)-unlist(lastp_cd))^2))

                  # Update parameter list.
                  lastp_cd <- p_cd

                  # Update the iteration number.
                  iter_cd = iter_cd + 1

                } # End coordinate descent looping.


              } # End Rasch conditional.

            } # End looping across covariates.

          } # End anchor item condtional.


        } # End item type conditional.

      } # End looping through items.


    } # End optimization method conditional.




    if(method == "MNR") {
      inv_hess_diag[[num_items+1]] <-
        inv_hess_impact_diag[1:ncol(mean_predictors)]
      inv_hess_diag[[num_items+2]] <-
        inv_hess_impact_diag[-(1:ncol(mean_predictors))]
    }



    if(method == "MNR" || method == "UNR") break


    p_cd_all <- p_cd

    # Update and check for convergence: Calculate the difference
    # in parameter estimates from current to previous.
    eps_cd_all = sqrt(sum((unlist(p_cd_all)-unlist(lastp_cd_all))^2))

    # Update parameter list.
    lastp_cd_all <- p_cd_all

    cat('\r', sprintf("CD Iteration: %d  CD Change: %f",
                      iter_cd_all,
                      round(eps_cd_all, nchar(final_control$tol))))
    utils::flush.console()

    # Update the iteration number.
    iter_cd_all = iter_cd_all + 1

    } # End CD iterations.


    if(max_tau) {

      id_max_z <- max(abs(id_max_z))
      return(id_max_z)

    } else {

      if(method == "MNR") {
        return(list(p=p,
                    inv_hess_diag=inv_hess_diag,
                    under_identified=under_identified))
      } else if(method == "UNR") {
        return(list(p=p,
                    under_identified=under_identified))
      } else if(method == "CD") {
        return(list(p=p_cd_all,
                    under_identified=under_identified))
      }

    }

  }
