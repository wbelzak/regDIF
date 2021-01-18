#' Identify the minimum value of tau such that all DIF is removed.
#'
#' @param p List of parameters with starting values obtained from preprocess.
#' @param item_data Matrix or dataframe of item responses.
#' @param pred_data Matrix or dataframe of DIF and/or impact predictors.
#' @param mean_predictors Possibly different matrix of predictors for the mean
#' impact equation.
#' @param var_predictors Possibly different matrix of predictors for the
#' variance impact equation.
#' @param item_type Optional character value or vector indicating the type of
#' item to be modeled.
#' @param theta Vector of fixed weights to approximate the latent variable.
#' @param pen_type Character value indicating the penalty function to use.
#' @param tau_vec Vector of tau values that either are automatically generated
#' or provided by the user.
#' @param num_tau Number of tau values to run Reg-DIF on.
#' @param alpha Numeric value indicating the alpha parameter in the elastic net
#' penalty function.
#' @param gamma Numeric value indicating the gamma parameter in the MCP
#' function.
#' @param anchor Optional numeric value or vector indicating which item
#' response(s) are anchors (e.g., \code{anchor = 1}).
#' @param final_control List of control parameters.
#' @param samp_size Numeric value indicating the sample size.
#' @param num_items Numeric value indicating the number of items.
#' @param num_responses Vector with number of responses for each item.
#' @param num_predictors Numeric value indicating the number of predictors.
#' @param num_quad Numeric value indicating the number of quadrature points.
#' @param adapt_quad Logical value indicating whether to use adaptive quad.
#' @param optim_method Character value indicating the type of optimization
#' method to use.
#' @param em_history List to save EM iterations for supplemental EM algorithm.
#'
#' @keywords internal
#'
identify_tau <- function(p,
                         item_data,
                         pred_data,
                         mean_predictors,
                         var_predictors,
                         item_type,
                         theta,
                         pen_type,
                         tau_vec,
                         num_tau,
                         alpha,
                         gamma,
                         anchor,
                         final_control,
                         samp_size,
                         num_items,
                         num_responses,
                         num_predictors,
                         num_quad,
                         adapt_quad,
                         optim_method,
                         em_history) {

  # Maximization settings.
  lastp <- p
  eps <- Inf
  iter <- 0

  # Loop until convergence or maximum number of iterations.
  while(eps > final_control$tol &&
        iter < final_control$maxit){

    # E-step: Evaluate Q function with current parameter estimates p.
    etable <- Estep(p,
                    item_data,
                    pred_data,
                    mean_predictors,
                    var_predictors,
                    theta,
                    samp_size,
                    num_items,
                    num_responses,
                    adapt_quad,
                    num_quad)

    if(optim_method == "multi") {
      # M-step: Optimize parameters using multivariate NR.
      p <- Mstep_block(p,
                       item_data,
                       pred_data,
                       mean_predictors,
                       var_predictors,
                       etable,
                       item_type,
                       pen_type,
                       tau_vec[1],
                       alpha,
                       gamma,
                       anchor,
                       samp_size,
                       num_responses,
                       num_items,
                       num_quad,
                       num_predictors)
    } else if(optim_method == "uni") {
      # M-step: Optimize parameters using one round of coordinate descent.
      p <- Mstep_cd(p,
                    item_data,
                    pred_data,
                    mean_predictors,
                    var_predictors,
                    etable,
                    item_type,
                    pen_type,
                    tau_vec[1],
                    alpha,
                    gamma,
                    anchor,
                    samp_size,
                    num_responses,
                    num_items,
                    num_quad,
                    num_predictors)
    }

    # Update and check for convergence: Calculate the difference
    # in parameter estimates from current to previous.
    eps = sqrt(sum((unlist(p)-unlist(lastp))^2))

    # Save parameter estimates for supplemental em.
    em_history[[1]][,iter] <- unlist(p)
    if(eps > final_control$tol) {
      em_history[[1]] <- cbind(em_history[[1]],
                               matrix(0,ncol=1,nrow=length(unlist(p))))
    }

    # Update parameter list.
    lastp <- p

    # Update the iteration number.
    iter = iter + 1
    if(iter == final_control$maxit) {
      warning(paste0("Could not find a suitable vector of tuning parameter",
                     "values. Please input your own vector of tau values."))
    }

    cat('\r', sprintf("Models Completed: %d of %d  Iteration: %d  Change: %f",
                      0,
                      num_tau,
                      iter,
                      round(eps, nchar(final_control$tol))))

    utils::flush.console()


  }

  # Get information criteria.
  infocrit <- information_criteria(etable,
                                   p,
                                   item_data,
                                   pred_data,
                                   mean_predictors,
                                   var_predictors,
                                   gamma,
                                   samp_size,
                                   num_responses,
                                   num_items,
                                   num_quad)
  # Final M-step.
  if(optim_method == "multi") {
    max_tau <- Mstep_block_idtau(p,
                                 item_data,
                                 pred_data,
                                 mean_predictors,
                                 var_predictors,
                                 etable,
                                 item_type,
                                 pen_type,
                                 tau_vec[1],
                                 alpha,
                                 gamma,
                                 anchor,
                                 samp_size,
                                 num_responses,
                                 num_items,
                                 num_quad,
                                 num_predictors)
  } else if(optim_method == "uni") {
    max_tau <- Mstep_cd_idtau(p,
                              item_data,
                              pred_data,
                              mean_predictors,
                              var_predictors,
                              etable,
                              item_type,
                              pen_type,
                              tau_vec[1],
                              alpha,
                              gamma,
                              anchor,
                              samp_size,
                              num_responses,
                              num_items,
                              num_quad,
                              num_predictors)
  }



  return(list(p=p,infocrit=infocrit,max_tau=max_tau,em_history=em_history))
}



#' Final M-step uing coordinate descent for identifying the minimum value
#' of tau.
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
#' @param num.quad Number of quadrature points used for approximating the
#' latent variable.
#' @param num_predictors Number of predictors.
#'
#' @keywords internal
#'
Mstep_cd_idtau <-
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

      if(any(item == anchor)) next

        p2 <- unlist(p)

        # Intercept DIF updates.
        for(cov in 1:num_predictors) {

          # End routine if only one anchor item is left on each covariate
          # for each item parameter.
          if(is.null(anchor) &&
             sum(p2[grep(paste0("c1(.*?)cov",cov),names(p2))] != 0) >
             (num_items - 2) &&
             alpha == 1){
            next
          }

          anl_deriv <- d_bernoulli("c1",
                                       p[[item]],
                                       etable_item,
                                       theta,
                                       pred_data,
                                       cov,
                                       samp_size,
                                       num_items,
                                       num_quad)
          z <- p[[item]][[2+cov]] - anl_deriv[[1]]/anl_deriv[[2]]
          id_max_z <- c(id_max_z,z)
        }

        if(item_type[item] != "Rasch") next

        # Slope DIF updates.
        for(cov in 0:(num_predictors-1)) {

          # End routine if only one anchor item is left on each covariate
          # for each item parameter.
          if(is.null(anchor) &
             sum(p2[grep(paste0("a1(.*?)cov",cov+1),names(p2))] != 0) >
             (num_items - 2) &
             alpha == 1){
            next
          }

            anl_deriv <- d_bernoulli("a1",
                                         p[[item]],
                                         etable_item,
                                         theta,
                                         pred_data,
                                         cov,
                                         samp_size,
                                         num_items,
                                         num_quad)
            z <- p[[item]][[2+num_predictors+cov]] -
              anl_deriv[[1]]/anl_deriv[[2]]
            id_max_z <- c(id_max_z,z)
          }




    } else if(num_responses[item] > 2) {

      if(!any(item == anchor)){

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
             (num_items - 2) &
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
