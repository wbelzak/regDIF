#' Maximization step.
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
Mstep <-
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


  # Obtain parameter estimates and posterior probabilities.
  p_impact <- c(p[[num_items+1]],p[[num_items+2]])
  etable_all <- elist[[2]]

  # Impact mean updates.
  mean_parms <- grep(paste0("g"),names(p_impact),fixed=T)
  var_parms <- grep(paste0("b"),names(p_impact),fixed=T)

  for(cov in 0:(ncol(mean_predictors)-1)) {
    anl_deriv <- d_alpha_cpp(p[[num_items+1]],
                             p[[num_items+2]],
                             etable_all,
                             elist$theta,
                             mean_predictors,
                             var_predictors,
                             cov=cov,
                             samp_size,
                             num_items,
                             num.quad)
    p_new <- p_impact[mean_parms][cov+1] - anl_deriv[[1]]/anl_deriv[[2]]
    p_impact <- replace(p_impact,names(p_new),p_new)
  }

  # Impact variance updates.
  for(cov in 0:(ncol(var_predictors)-1)) {
    anl_deriv <- d_phi_cpp(p[[num_items+1]],
                           p[[num_items+2]],
                           etable_all,
                           elist$theta,
                           mean_predictors,
                           var_predictors,
                           cov=cov,
                           samp_size,
                           num_items,
                           num.quad)
    p_new <- p_impact[var_parms][cov+1] - anl_deriv[[1]]/anl_deriv[[2]]
    p_impact <- replace(p_impact,names(p_new),p_new)
  }

  g <- p_impact[mean_parms]
  b <- p_impact[var_parms]

  p[[num_items+1]] <- replace(p[[num_items+1]],names(g),g)
  p[[num_items+2]] <- replace(p[[num_items+2]],names(b),b)

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

      # Intercept updates.
      c0_parms <- grep(paste0("c0_itm",item,"_"),names(p_item),fixed=T)
      anl_deriv <- d_bernoulli_cpp("c0",
                                   p_item,
                                   etable[[1]],
                                   etable[[2]],
                                   elist$theta,
                                   pred.data,
                                   cov=0,
                                   samp_size,
                                   num_items,
                                   num.quad)
      p_new <- p_item[c0_parms][1] - anl_deriv[[1]]/anl_deriv[[2]]
      p_item <- replace(p_item,names(p_new),p_new)

      # Slope updates.
      if(item.type[item] != "Rasch") {
        a0_parms <- grep(paste0("a0_itm",item,"_"),names(p_item),fixed=T)
        anl_deriv <- d_bernoulli_cpp("a0",
                                     p_item,
                                     etable[[1]],
                                     etable[[2]],
                                     elist$theta,
                                     pred.data,
                                     cov=0,
                                     samp_size,
                                     num_items,
                                     num.quad)
        p_new <- p_item[a0_parms] - anl_deriv[[1]]/anl_deriv[[2]]
        p_item <- replace(p_item,names(p_new),p_new)
      }


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
          p_new <- ifelse(pen.type == "lasso",
                          soft_thresh_cpp(z,alpha,tau_current),
                          firm_thresh_cpp(z,alpha,tau_current,gamma))
          names(p_new) <- names(z)
          p_item <- replace(p_item,names(p_new),p_new)
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
            p_new <- ifelse(pen.type == "mcp",
                            firm_thresh_cpp(z,alpha,tau_current,gamma),
                            soft_thresh_cpp(z,alpha,tau_current))
            names(p_new) <- names(z)
            p_item <- replace(p_item,names(p_new),p_new)
          }

        }

      }

    } else if(num_responses[item] > 2) {

      c0_parms <- grep(paste0("c0_itm",item,"_"),names(p_item),fixed=T)

      # Intercept updates.
      anl_deriv <- d_categorical("c0",
                                     p_item,
                                     etable,
                                     elist$theta,
                                     pred.data,
                                     thr=-1,
                                     cov=-1,
                                     samp_size,
                                     num_responses[[item]],
                                     num_items,
                                 num.quad)
      p_new <- p_item[c0_parms][1] - anl_deriv[[1]]/anl_deriv[[2]]
      p_item <- replace(p_item,names(p_new),p_new)

      # Threshold updates.
      if(num_responses[item] > 2) {
        for(thr in 2:(num_responses[item]-1)) {
          anl_deriv <- d_categorical("c0",
                                         p_item,
                                         etable,
                                         elist$theta,
                                         pred.data,
                                         thr=thr,
                                         cov=-1,
                                         samp_size,
                                         num_responses[[item]],
                                         num_items,
                                     num.quad)
          p_new <- p_item[c0_parms][thr] - anl_deriv[[1]]/anl_deriv[[2]]
          p_item <- replace(p_item,names(p_new),p_new)
        }
      }

      # Slope updates.
      if(item.type[item] != "Rasch") {
        a0_parms <- grep(paste0("a0_itm",item,"_"),names(p_item),fixed=T)
        anl_deriv <- d_categorical("a0",
                                       p_item,
                                       etable,
                                       elist$theta,
                                       pred.data,
                                       thr=-1,
                                       cov=-1,
                                       samp_size,
                                       num_responses[[item]],
                                       num_items,
                                   num.quad)
        p_new <- p_item[a0_parms] - anl_deriv[[1]]/anl_deriv[[2]]
        p_item <- replace(p_item,names(p_new),p_new)
      }

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
          p_new <- ifelse(pen.type == "mcp",
                          firm_thresh_cpp(z,alpha,tau_current,gamma),
                          soft_thresh_cpp(z,alpha,tau_current))
          names(p_new) <- names(z)
          p_item <- replace(p_item,names(p_new),p_new)
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
            p_new <- ifelse(pen.type == "mcp",
                            firm_thresh_cpp(z,alpha,tau_current,gamma),
                            soft_thresh_cpp(z,alpha,tau_current))
            names(p_new) <- names(z)
            p_item <- replace(p_item,names(p_new),p_new)
          }
        }
      }


      # Gaussian responses.
    } else if(num_responses[item] == 1) {

      # Intercept updates.
      c0_parms <- grep(paste0("c0_itm",item,"_"),names(p_item),fixed=T)
      anl_deriv <- d_mu_gaussian("c0",
                                 p_item,
                                 etable,
                                 elist$theta,
                                 item.data[,item],
                                 pred.data,
                                 cov=NULL,
                                 samp_size,
                                 num_items,
                                 num.quad)
      p_new <- p_item[c0_parms][1] - anl_deriv[[1]]/anl_deriv[[2]]
      p_item <- replace(p_item,names(p_new),p_new)

      # Slope updates.
      if(item.type[item] != "Rasch") {
        a0_parms <- grep(paste0("a0_itm",item,"_"),names(p_item),fixed=T)
        anl_deriv <- d_mu_gaussian("a0",
                                   p_item,
                                   etable,
                                   elist$theta,
                                   item.data[,item],
                                   pred.data,
                                   cov=NULL,
                                   samp_size,
                                   num_items,
                                   num.quad)
        p_new <- p_item[a0_parms] - anl_deriv[[1]]/anl_deriv[[2]]
        p_item <- replace(p_item,names(p_new),p_new)
      }

      # Residual updates.
      s0_parms <- grep(paste0("s0_itm",item,"_"),names(p_item),fixed=T)
      anl_deriv <- d_sigma_gaussian("s0",
                                    p_item,
                                    etable,
                                    elist$theta,
                                    item.data[,item],
                                    pred.data,
                                    cov=NULL,
                                    samp_size,
                                    num_items,
                                    num.quad)
      p_new <- p_item[s0_parms][1] - anl_deriv[[1]]/anl_deriv[[2]]
      p_item <- replace(p_item,names(p_new),p_new)


      if(!any(item == anchor)) {

        # Residual DIF updates.
        for(cov in 1:num_predictors) {
          s1_parms <-
            grep(paste0("s1_itm",item,"_cov",cov),names(p_item),fixed=T)
          anl_deriv <- d_sigma_gaussian("s1",
                                        p_item,
                                        etable,
                                        elist$theta,
                                        item.data[,item],
                                        pred.data,
                                        cov=cov,
                                        samp_size,
                                        num_items,
                                        num.quad)
          p_new <- p_item[s1_parms][1] - anl_deriv[[1]]/anl_deriv[[2]]
          p_item <- replace(p_item,names(p_new),p_new)
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
          p_new <- ifelse(pen.type == "mcp",
                          firm_thresh_cpp(z,alpha,tau_current,gamma),
                          soft_thresh_cpp(z,alpha,tau_current))
          names(p_new) <- names(z)
          p_item <- replace(p_item,names(p_new),p_new)
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
            p_new <- ifelse(pen.type == "mcp",
                            firm_thresh_cpp(z,alpha,tau_current,gamma),
                            soft_thresh_cpp(z,alpha,tau_current))
            names(p_new) <- names(z)
            p_item <- replace(p_item,names(p_new),p_new)
          }
        }
      }

    }


    p[[item]] <- replace(p[[item]],names(p_item),p_item)


  }

  return(p)

}

