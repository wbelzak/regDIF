#' Identify the minimum value of tau such that all DIF is removed.
#'
#' @param p List of parameters with starting values obtained from preprocess.
#' @param item.data Matrix or dataframe of item responses.
#' @param pred.data Matrix or dataframe of DIF and/or impact predictors.
#' @param mean_predictors Possibly different matrix of predictors for the mean
#' impact equation.
#' @param var_predictors Possibly different matrix of predictors for the
#' variance impact equation.
#' @param item.type Optional character value or vector indicating the type of
#' item to be modeled.
#' @param pen.type Character value indicating the penalty function to use.
#' @param tau_vec Vector of tau values that either are automatically generated
#' or provided by the user.
#' @param alpha Numeric value indicating the alpha parameter in the elastic net
#' penalty function.
#' @param gamma Numeric value indicating the gamma parameter in the MCP
#' function.
#' @param anchor Optional numeric value or vector indicating which item
#' response(s) are anchors (e.g., \code{anchor = 1}).
#' @param final.control List of control parameters.
#' @param samp_size Numeric value indicating the sample size.
#' @param num_items Numeric value indicating the number of items.
#' @param num_responses Vector with number of responses for each item.
#' @param num_predictors Numeric value indicating the number of predictors.
#' @param num.quad Numeric value indicating the number of quadrature points.
#' @param adapt.quad Logical value indicating whether to use adaptive quad.
#'
#' @keywords internal
#'
identify_tau <- function(p,
                         item.data,
                         pred.data,
                         mean_predictors,
                         var_predictors,
                         item.type,
                         pen.type,
                         tau_vec,
                         alpha,
                         gamma,
                         anchor,
                         final.control,
                         samp_size,
                         num_items,
                         num_responses,
                         num_predictors,
                         num.quad,
                         adapt.quad) {

  # Maximization settings.
  lastp <- p
  eps <- Inf
  iter <- 0

  # Loop until convergence or maximum number of iterations.
  while(eps > final.control$tol &&
        iter < final.control$maxit){

    # E-step: Evaluate Q function with current parameter estimates p.
    elist <- Estep(p,
                   item.data,
                   pred.data,
                   mean_predictors,
                   var_predictors,
                   samp_size,
                   num_items,
                   num_responses,
                   adapt.quad,
                   num.quad)

    # M-step: Optimize parameters.
    p <- Mstep(p,
               item.data,
               pred.data,
               mean_predictors,
               var_predictors,
               elist,
               item.type,
               pen.type,
               tau_vec[1],
               alpha,
               gamma,
               anchor,
               samp_size,
               num_responses,
               num_items,
               num.quad,
               num_predictors)

    # Update and check for convergence: Calculate the difference
    # in parameter estimates from current to previous.
    eps = sqrt(sum((unlist(p)-unlist(lastp))^2))

    # Update parameter list.
    lastp <- p

    # Update the iteration number.
    iter = iter + 1
    if(iter == final.control$maxit) {
      warning(paste0("Could not find a suitable vector of tuning parameter",
                     "values. Please input your own vector of tau values."))
    }

    cat('\r', sprintf("Generating vector of tuning parameters..."))

    utils::flush.console()


  }

  max_tau <- Mstep_id_tau(p,
                          item.data,
                          pred.data,
                          mean_predictors,
                          var_predictors,
                          elist,
                          item.type,
                          pen.type,
                          tau_vec[1],
                          alpha,
                          gamma,
                          anchor,
                          samp_size,
                          num_responses,
                          num_items,
                          num.quad,
                          num_predictors)

  return(list(p=p,max_tau=max_tau))
}



#' Final M-step for identifying the minimum value of tau.
#'
#' @param p List of parameters.
#' @param item.data Matrix or data frame of item responses.
#' @param pred.data Matrix or data frame of DIF and/or impact predictors.
#' @param mean_predictors Possibly different matrix of predictors for the mean
#' impact equation.
#' @param var_predictors Possibly different matrix of predictors for the
#' variance impact equation.
#' @param elist List of E-tables for item and impact equations, in addition to
#' theta values.
#' @param item.type Optional character value or vector indicating the type of
#' item to be modeled.
#' @param pen.type Character value indicating the penalty function to use.
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
Mstep_id_tau <-
  function(p,
           item.data,
           pred.data,
           mean_predictors,
           var_predictors,
           elist,
           item.type,
           pen.type,
           tau_current,
           alpha,
           gamma,
           anchor,
           samp_size,
           num_responses,
           num_items,
           num.quad,
           num_predictors) {


  # Last Mstep
  id_max_z <- 0

  # Maximize independent logistic regressions.
  for (item in 1:num_items) {

    # Get posterior probabilities.
    etable <- replicate(n=num_responses[item],
                        elist[[1]][[item]],
                        simplify = F)

    # Obtain E-tables for each response category.
    if(num_responses[item] > 1) {
      for(resp in 1:num_responses[item]) {
        etable[[resp]][which(
          !(etable[[resp]][, ncol(etable[[resp]])] == resp)), ] <- 0
      }
    }

    # Get item parameters.
    p_item <- p[[item]]
    etable <- lapply(etable, function(x) x[,1:num.quad])

    # Bernoulli responses.
    if(num_responses[item] == 2) {

      if(!any(item == anchor)) {

        p2 <- unlist(p)

        # Intercept DIF updates.
        for(cov in 0:(num_predictors-1)) {

          # End routine if only one anchor item is left on each covariate
          # for each item parameter.
          if(is.null(anchor) &&
             sum(p2[grep(paste0("c1(.*?)cov",cov+1),names(p2))] != 0) >
             (num_items - 2) &&
             alpha == 1){
            next
          }

          c1_parms <-
            grep(paste0("c1_itm",item,"_cov",cov+1),names(p_item),fixed=T)
          anl_deriv <- d_bernoulli_cpp("c1",
                                       p_item,
                                       etable[[1]],
                                       etable[[2]],
                                       elist$theta,
                                       pred.data,
                                       cov,
                                       samp_size,
                                       num_items,
                                       num.quad)
          z <- p_item[c1_parms] - anl_deriv[[1]]/anl_deriv[[2]]
          id_max_z <- c(id_max_z,z)
        }

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

          if(item.type[item] != "Rasch") {
            a1_parms <-
              grep(paste0("a1_itm",item,"_cov",cov+1),names(p_item),fixed=T)
            anl_deriv <- d_bernoulli_cpp("a1",
                                         p_item,
                                         etable[[1]],
                                         etable[[2]],
                                         elist$theta,
                                         pred.data,
                                         cov,
                                         samp_size,
                                         num_items,
                                         num.quad)
            z <- p_item[a1_parms] - anl_deriv[[1]]/anl_deriv[[2]]
            id_max_z <- c(id_max_z,z)
          }

        }

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

          c1_parms <-
            grep(paste0("c1_itm",item,"_cov",cov),names(p_item),fixed=T)
          anl_deriv <- d_categorical("c1",
                                     p_item,
                                     etable,
                                     elist$theta,
                                     pred.data,
                                     thr=-1,
                                     cov,
                                     samp_size,
                                     num_responses[[item]],
                                     num_items,
                                     num.quad)
          z <- p_item[c1_parms] - anl_deriv[[1]]/anl_deriv[[2]]
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

          if(item.type[item] != "Rasch") {
            a1_parms <-
              grep(paste0("a1_itm",item,"_cov",cov),names(p_item),fixed=T)
            anl_deriv <- d_categorical("a1",
                                       p_item,
                                       etable,
                                       elist$theta,
                                       pred.data,
                                       thr=-1,
                                       cov,
                                       samp_size,
                                       num_responses[[item]],
                                       num_items,
                                       num.quad)
            z <- p_item[a1_parms] - anl_deriv[[1]]/anl_deriv[[2]]
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
            grep(paste0("c1_itm",item,"_cov",cov),names(p_item),fixed=T)
          anl_deriv <- d_mu_gaussian("c1",
                                     p_item,
                                     etable,
                                     elist$theta,
                                     item.data[,item],
                                     pred.data,
                                     cov,
                                     samp_size,
                                     num_items,
                                     num.quad)
          z <- p_item[c1_parms] - anl_deriv[[1]]/anl_deriv[[2]]
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

          if(item.type[item] != "Rasch"){
            a1_parms <-
              grep(paste0("a1_itm",item,"_cov",cov),names(p_item),fixed=T)
            anl_deriv <- d_mu_gaussian("a1",
                                       p_item,
                                       etable,
                                       elist$theta,
                                       item.data[,item],
                                       pred.data,
                                       cov,
                                       samp_size,
                                       num_items,
                                       num.quad)
            z <- p_item[a1_parms] - anl_deriv[[1]]/anl_deriv[[2]]
            id_max_z <- c(id_max_z,z)
          }
        }
      }

    }



  }


  id_max_z <- max(abs(id_max_z))
  return(id_max_z)

}
