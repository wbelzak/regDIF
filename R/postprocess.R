#' Maximization step.
#'
#' @param estimates List of converged parameters.
#' @param item.data User-given matrix or data.frame of DIF and/or impact
#' predictors.
#' @param pred.data User-given matrix or data.frame of item responses.
#' @param item_data Processed matrix or data.frame of item responses.
#' @param pred_data Processed matrix or data.frame of DIF and/or impact
#' predictors.
#' @param mean_predictors Possibly different matrix of predictors for the mean
#' impact equation.
#' @param var_predictors Possibly different matrix of predictors for the
#' variance impact equation.
#' @param tau_vec Optional numeric vector of tau values.
#' @param num_tau Logical indicating whether the minimum tau value needs to be
#' identified during the regDIF procedure.
#' @param alpha Numeric value indicating the alpha parameter in the elastic net
#' penalty function.
#' @param pen Tuning parameter index.
#' @param anchor Anchor item(s).
#' @param control Optional list of user-defined control parameters
#' @param final_control List of final control parameters.
#' @param final List of model results.
#' @param samp_size Sample size in dataset.
#' @param num_responses Number of responses for each item.
#' @param num_predictors Number of predictors.
#' @param num_items Number of items in dataset.
#' @param num_quad Number of quadrature points used for approximating the
#' latent variable.
#'
#' @keywords internal
#'
postprocess <-
  function(estimates,
           item.data,
           pred.data,
           item_data,
           pred_data,
           mean_predictors,
           var_predictors,
           tau_vec,
           num_tau,
           alpha,
           pen,
           anchor,
           control,
           final_control,
           final,
           samp_size,
           num_responses,
           num_predictors,
           num_items,
           num_quad) {

  # Get estimates and information criteria.
  p <- estimates$p
  infocrit <- estimates$infocrit
  em_history <- estimates$em_history
  complete_info <- estimates$complete_info
  under_identified <- estimates$under_identified

  # Organize impact parameters.
  if(is.null(control$impact.mean.data)) {
    if(is.null(colnames(pred_data)) |
       length(colnames(pred_data)) == 0) {
      mean_names <- paste0('mean.cov',1:ncol(mean_predictors))
    } else {
      mean_names <- paste0('mean.',colnames(pred_data))
    }

  } else {
    if(is.null(colnames(control$impact.mean.data)) |
       length(colnames(control$impact.mean.data)) == 0) {
      mean_names <- paste0('mean.cov',1:ncol(mean_predictors))
    } else {
      mean_names <- paste0('mean.',colnames(control$impact.mean.data))
    }

  }

  if(is.null(control$impact.var.data)) {
    if(is.null(colnames(pred_data)) |
       length(colnames(pred_data)) == 0) {
      var_names <- paste0('var.cov',1:ncol(var_predictors))
    } else {
      var_names <- paste0('var.',colnames(pred_data))
    }

  } else {
    if(is.null(colnames(pred_data)) |
       length(colnames(pred_data)) == 0){
      var_names <- paste0('var.cov',1:ncol(var_predictors))
    } else {
      var_names <- paste0('var.',colnames(control$impact.var.data))
    }

  }

  lv_parms <- c(p[[num_items+1]],p[[num_items+2]])
  lv_names <- c(mean_names,var_names)

  # Organize item baseline parameters.
  p2 <- unlist(p)
  all_items_parms_base <- NULL
  all_items_names_base <- NULL
  if(is.null(colnames(item_data)) |
     length(colnames(item_data)) == 0){
    item_names <- paste0("item",1:num_items)

  } else{
    item_names <- colnames(item_data)

  }

  for(item in 1:num_items) {
    if(num_responses[item] == 1) {
      item_parms_base <- c(p2[grep(paste0("c0_item",item,"_"),names(p2))],
                           p2[grep(paste0("a0_item",item,"_"),names(p2))],
                           p2[grep(paste0("s0_item",item,"_"),names(p2))])
      item_names_base <- c(paste0(item_names[item],".int"),
                           paste0(item_names[item],".slp"),
                           paste0(item_names[item],".res"))

    } else if(num_responses[item] == 2) {
      item_parms_base <- c(p2[grep(paste0("c0_item",item,"_"),names(p2))],
                           p2[grep(paste0("a0_item",item,"_"),names(p2))])
      item_names_base <- c(paste0(item_names[item],".int"),
                           paste0(item_names[item],".slp"))

    } else if(num_responses[item] > 2) {
      item_parms_base <- c(p2[grep(paste0("c0_item",item,"_"),names(p2))],
                           p2[grep(paste0("a0_item",item,"_"),names(p2))])
      item_names_base <- c(paste0(item_names[item],".int",
                                  1:(num_responses[item]-1)),
                           paste0(item_names[item],".slp"))

    }
    all_items_parms_base <- c(all_items_parms_base,item_parms_base)
    all_items_names_base <- c(all_items_names_base,item_names_base)

  }

  # Transform threshold values if ordered categorical item
  for(item in 1:num_items) {
    if(num_responses[item] > 2) {
      threshold_parms <- p2[grep(paste0("c0_item",item,"_"),
                            names(p2))][2:(num_responses[item]-1)]
      intercept_parm <- p2[grep(paste0("c0_item",item,"_"),names(p2))][1]

      all_items_parms_base[grep(paste0("c0_item",item,"_"),
                                names(all_items_parms_base))][2:(
                                  num_responses[item]-1)] <-
        intercept_parm - threshold_parms
    } else {
      next
    }
  }



  # Organize item DIF parameters.
  all_items_parms_dif <- NULL
  all_items_names_dif <- NULL
  if(is.null(colnames(pred_data)) |
     length(colnames(pred_data)) == 0) {
    cov_names <- paste0("cov",1:num_predictors)

  } else {
    cov_names <- colnames(pred_data)

  }

  for(item in 1:num_items) {
    if(num_responses[item] == 1) {
      item_parms_dif <- c(p2[grep(paste0("c1_item",item,"_"),names(p2),fixed=T)],
                          p2[grep(paste0("a1_item",item,"_"),names(p2),fixed=T)],
                          p2[grep(paste0("s1_item",item,"_"),names(p2),fixed=T)])
      item_names_dif <- c(paste0(rep(item_names[item],
                                     each = num_predictors),'.int.',cov_names),
                          paste0(rep(item_names[item],
                                     each = num_predictors),'.slp.',cov_names),
                          paste0(rep(item_names[item],
                                     each = num_predictors),'.res.',cov_names))

    } else {
      item_parms_dif <- c(p2[grep(paste0("c1_item",item,"_"),names(p2),fixed=T)],
                          p2[grep(paste0("a1_item",item,"_"),
                                  names(p2),
                                  fixed=T)])
      item_names_dif <- c(paste0(rep(item_names[item],
                                     each = num_predictors),'.int.',cov_names),
                          paste0(rep(item_names[item],
                                     each = num_predictors),'.slp.',cov_names))

    }
    all_items_parms_dif <- c(all_items_parms_dif,item_parms_dif)
    all_items_names_dif <- c(all_items_names_dif,item_names_dif)

  }

  # If model becomes under-identified, return NA values.
  if(under_identified) {
    # Print information about optimization.
    cat('\r',
        sprintf(paste0("Models Completed: %d of %d  Iteration: %d  Change: %d",
                       "              "),
                pen, length(tau_vec), 0, 0))

    if(pen != length(tau_vec)) utils::flush.console()
    return(final)
    }

  # Assign rest of output to final list.
  final$tau_vec[pen] <- tau_vec[pen]
  final$aic[pen] <- round(infocrit$aic,3)
  final$bic[pen] <- round(infocrit$bic,3)
  final$impact[,pen] <- round(lv_parms,3)
  final$base[,pen] <- round(all_items_parms_base,3)
  final$dif[,pen] <- round(all_items_parms_dif,3)
  final$em_history[[pen]] <- em_history[[pen]]
  final$complete_ll_info <- complete_info
  # final$log_like[,pen] <- c(round(infocrit$complete_ll,3),
  #                           round(infocrit$observed_ll,3))
  final$data <- list(item.data=item.data, pred.data=pred.data)
  rownames(final$impact) <- lv_names
  rownames(final$base) <- all_items_names_base
  rownames(final$dif) <- all_items_names_dif
  if(!(any(num_responses > 2))) {
    for(item in 1:num_items) {
      names(final$complete_ll_info[[item]]) <- names(p[[item]])
    }
    names(final$complete_ll_info[[num_items+1]]) <- names(p[[num_items+1]])
    names(final$complete_ll_info[[num_items+2]]) <- names(p[[num_items+2]])
  }
  # rownames(final$log_like) <- c("complete","observed")

  # Order item parms.
  final_int_thr_base <-
    final$base[grep(paste0(c(".int",".thr"),
                           collapse = "|"),
                    rownames(final$base)),
               pen]
  final_slp_base <-
    final$base[grep(".slp",
                    rownames(final$base)),
               pen]
  final_res_base <-
    final$base[grep(".res",
                    rownames(final$base)),
               pen]
  final_names_base <-
    names(c(final_int_thr_base,final_slp_base,final_res_base))
  final$base[,pen] <-
    matrix(c(final_int_thr_base,final_slp_base,final_res_base), ncol = 1)
  rownames(final$base) <- final_names_base

  final_int_dif <- final$dif[grep(".int",
                                  rownames(final$dif)),
                             pen]
  final_slp_dif <- final$dif[grep(".slp",
                                  rownames(final$dif)),
                             pen]
  final_res_dif <- final$dif[grep(".res",
                                  rownames(final$dif)),
                             pen]
  final_names_dif <- names(c(final_int_dif,final_slp_dif,final_res_dif))
  final$dif[,pen] <-
    matrix(c(final_int_dif,final_slp_dif,final_res_dif), ncol = 1)
  rownames(final$dif) <- final_names_dif

  # Warn if there is a large change in DIF parameters.
  if(pen > 1){
    second_last <- sum(final$dif[-1,pen-1] == 0)
    last <- sum(final$dif[-1,pen] == 0)

    if((second_last - last) > (num_predictors*num_items)) {
      warning(paste0("Large increase in the number of DIF parameters ",
                  "from iteration ",
                  pen-1,
                  " to ",
                  pen,
                  ".\n  Two Options:\n  1. Provide smaller differences ",
                  "between tau values.\n  2. Provide anchor item(s)."),
           call. = FALSE)
    }

  }

  # Get parameter estimates to determine if tau value is too small to remove all
  # dif from model.
  p2 <- unlist(estimates$p)
  dif_parms <- p2[grep(paste0("cov"),names(p2))]

  # Warn if tau does not remove all DIF
  if(is.null(anchor) &&
     pen == 1 &&
     sum(abs(dif_parms)) > 0 &&
     alpha == 1 &&
     num_tau >= 10) {
    warning(paste0("\nAutomatically-generated or user-defined ",
                   "tau value is too small to penalize all parameters to ",
                   "zero without anchor item. Larger values of tau are ",
                   "recommended."))
  }

  # Print information about optimization.
  cat('\r',
      sprintf(paste0("Models Completed: %d of %d  Iteration: %d  Change: %d",
                     "              "),
              pen, length(tau_vec), 0, 0))

  if(pen != length(tau_vec)) utils::flush.console()

  return(final)

}
