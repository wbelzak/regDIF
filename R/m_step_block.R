#' Maximization step using latent variable and item response blocks.
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
#'
#' @return a \code{"list"} of estimates obtained from the maximization step using multivariate
#' Newton-Raphson
#'
#' @importFrom foreach %dopar%
#'
#' @keywords internal
#'
Mstep_block <-
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
           max_tau) {

    # Set under-identified model to FALSE until proven TRUE.
    under_identified <- FALSE

    # Last Mstep
    if(max_tau) id_max_z <- 0

    # Latent impact updates.
    if(is.null(prox_data)) {
      anl_deriv_impact <- d_impact_block(p[[num_items+1]],
                                         p[[num_items+2]],
                                         eout$etable,
                                         eout$theta,
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

    # Item response updates.
    if(final_control$parallel[[1]]) {
      parallel::clusterExport(final_control$parallel[[2]], c("num_responses", "pred_data", "item_data",
                          "samp_size", "num_predictors", "num_quad", "num_items",
                          "prox_data", "item_type", "anchor", "pen_type",
                          "num_tau", "pen", "tau_vec", "alpha", "final_control",
                          "max_tau", "d_bernoulli_itemblock", "d_bernoulli_itemblock_proxy",
                          "grp_soft_threshold", "grp_firm_threshold", "soft_threshold",
                          "firm_threshold", "p", "eout", "inv_hess_diag", "tau_current",
                          "under_identified", "id_max_z"),
                    envir=environment())
      p_items <- foreach::foreach(item=1:num_items) %dopar% {

        if(num_responses[item] == 2) {

          if(is.null(prox_data)) {
            anl_deriv_item <- d_bernoulli_itemblock(p[[item]],
                                                    eout$etable,
                                                    eout$theta,
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
          if(any(item == anchor)) {
            return(list(unlist(p[item]),id_max_z))
          }

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
                 num_tau >= 10) {
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

            }

            # End routine if only one anchor item is left on each covariate
            # for each item parameter.
            if(is.null(anchor) &&
               sum(p2[grep(paste0("c1(.*?)cov",cov),names(p2))] != 0) >
               (num_items - 1) &&
               alpha == 1 &&
               (length(final_control$start.values) == 0 || pen > 1) &&
               num_tau >= 10) {
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
            if(item_type[item] == "rasch") {
              return(list(unlist(p[item]),id_max_z))
            }

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

          }
          if(under_identified) return(list(unlist(p[item]),id_max_z))


        }
        return(list(unlist(p[item]),id_max_z))

    }

      for(item in 1:num_items) {
        p[[item]] <- p_items[[item]][[1]]
      }

      } else{

        for (item in 1:num_items) {

          # Bernoulli responses.
          if(num_responses[item] == 2) {

            if(is.null(prox_data)) {
              anl_deriv_item <- d_bernoulli_itemblock(p[[item]],
                                                      eout$etable,
                                                      eout$theta,
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

              }

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

            }
            if(under_identified) break

          }
          }

        }





    inv_hess_diag[[num_items+1]] <-
      inv_hess_impact_diag[1:ncol(mean_predictors)]
    inv_hess_diag[[num_items+2]] <-
      inv_hess_impact_diag[-(1:ncol(mean_predictors))]

    if(max_tau) {
      if(final_control$parallel[[1]]) {
        id_max_z <- max(abs(unlist(sapply(1:num_items, function(items) p_items[[items]][[2]]))))
      } else{
        id_max_z <- max(abs(id_max_z))
      }
      return(id_max_z)
    } else {
      return(list(p=p,
                  inv_hess_diag=inv_hess_diag,
                  under_identified=under_identified))
    }

  }
