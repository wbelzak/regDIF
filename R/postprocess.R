###############
# Postprocess #
###############

postprocess <-
  function(estimates,
           responses,
           predictors,
           mean_predictors,
           var_predictors,
           y,
           x,
           impact.x,
           theta,
           lambda,
           alpha,
           pen,
           anchor,
           final.control,
           final,
           samp_size,
           num_responses,
           num_predictors,
           num_items,
           num_quadpts) {

    # responses <- data_scrub$responses;predictors <- data_scrub$predictors;mean_predictors <- data_scrub$mean_predictors;var_predictors <- data_scrub$var_predictors;theta <- data_scrub$theta;itemtypes <- data_scrub$itemtypes;lambda <- data_scrub$lambda; final.control <- data_scrub$final.control;final <- data_scrub$final;samp_size <- data_scrub$samp_size;num_items <- data_scrub$num_items;num_responses <- data_scrub$num_responses;num_predictors <- data_scrub$num_predictors;num_quadpts <- data_scrub$num_quadpts

  #get estimates and information criteria
  p <- estimates[[1]]
  infocrit <- estimates[[2]]

  #Organize impact parameters
  if(is.null(impact.x$mean)){ #mean
    if(is.null(colnames(x)) | length(colnames(x)) == 0){
      mean_names <- paste0('mean.cov',1:ncol(mean_predictors))
    } else{
      mean_names <- paste0('mean.',colnames(x))
    }
  } else{
    if(is.null(colnames(impact.x$mean)) | length(colnames(impact.x$mean)) == 0){
      mean_names <- paste0('mean.cov',1:ncol(mean_predictors))
    } else{
      mean_names <- paste0('mean.',colnames(impact.x$mean))
    }
  }
  if(is.null(impact.x$var)){
    if(is.null(colnames(x)) | length(colnames(x)) == 0){
      var_names <- paste0('var.cov',1:ncol(var_predictors))
    } else{
      var_names <- paste0('var.',colnames(x))
    }
  } else{
    if(is.null(colnames(x)) | length(colnames(x)) == 0){
      var_names <- paste0('var.cov',1:ncol(var_predictors))
    } else{
      var_names <- paste0('var.',colnames(impact.x$var))
    }
  }

  lv_parms <- c(p[[num_items+1]],p[[num_items+2]])
  lv_names <- c(mean_names,var_names)
  # if(is.null(colnames(x)) | length(colnames(x)) == 0){
  #   lv_names <- c(paste0('mean.cov',1:ncol(mean_predictors)),paste0('var.cov',1:ncol(var_predictors)))
  # } else{
  #   lv_names <- c(paste0('mean.',colnames(x)),paste0('var.',colnames(x)))
  # }

  #Organize item baseline parameters
  p2 <- unlist(p)
  all_items_parms_base <- NULL
  all_items_names_base <- NULL
  if(is.null(colnames(y)) | length(colnames(y)) == 0){
    item_names <- paste0("item",1:num_items)
  } else{
    item_names <- colnames(y)
  }

  for(item in 1:num_items){
    if(num_responses[item] == 1){
      item_parms_base <- c(p2[grep(paste0("c0_itm",item),names(p2))],p2[grep(paste0("a0_itm",item),names(p2))],p2[grep(paste0("s0_itm",item),names(p2))])
      item_names_base <- c(paste0(item_names[item],".int"),paste0(item_names[item],".slp"),paste0(item_names[item],".res"))
    } else if(num_responses[item] == 2){
      item_parms_base <- c(p2[grep(paste0("c0_itm",item),names(p2))],p2[grep(paste0("a0_itm",item),names(p2))])
      item_names_base <- c(paste0(item_names[item],".int"),paste0(item_names[item],".slp"))
    } else if(num_responses[item] > 2){
      item_parms_base <- c(p2[grep(paste0("c0_itm",item),names(p2))],p2[grep(paste0("a0_itm",item),names(p2))])
      item_names_base <- c(paste0(item_names[item],".int"),paste0(item_names[item],".thr",1:(num_responses[item]-2)),paste0(item_names[item],".slp"))
    }
    all_items_parms_base <- c(all_items_parms_base,item_parms_base)
    all_items_names_base <- c(all_items_names_base,item_names_base)
  }



  #Organize item dif parameters
  all_items_parms_dif <- NULL
  all_items_names_dif <- NULL
  if(is.null(colnames(x)) | length(colnames(x)) == 0){
    cov_names <- paste0("cov",1:num_predictors)
  } else{
    cov_names <- colnames(x)
  }

  for(item in 1:num_items){
    if(num_responses[item] == 1){
      item_parms_dif <- c(p2[grep(paste0("c1_itm",item),names(p2),fixed=T)],p2[grep(paste0("a1_itm",item),names(p2),fixed=T)],p2[grep(paste0("s1_itm",item),names(p2),fixed=T)])
      item_names_dif <- c(paste0(rep(item_names[item],each = num_predictors),'.int.',cov_names),paste0(rep(item_names[item],each = num_predictors),'.slp.',cov_names),paste0(rep(item_names[item],each = num_predictors),'.res.',cov_names))
    } else{
      item_parms_dif <- c(p2[grep(paste0("c1_itm",item),names(p2),fixed=T)],p2[grep(paste0("a1_itm",item),names(p2),fixed=T)])
      item_names_dif <- c(paste0(rep(item_names[item],each = num_predictors),'.int.',cov_names),paste0(rep(item_names[item],each = num_predictors),'.slp.',cov_names))
    }
    all_items_parms_dif <- c(all_items_parms_dif,item_parms_dif)
    all_items_names_dif <- c(all_items_names_dif,item_names_dif)
  }


  #assign output to final list
  final$lambda[pen] <- lambda[pen]
  final$aic[pen] <- round(infocrit[1],3)
  final$bic[pen] <- round(infocrit[2],3)
  final$impact.lv.parms[,pen] <- round(lv_parms,3)
  final$base.item.parms[,pen] <- round(all_items_parms_base,3)
  final$dif.item.parms[,pen] <- round(all_items_parms_dif,3)
  rownames(final$impact.lv.parms) <- lv_names
  rownames(final$base.item.parms) <- all_items_names_base
  rownames(final$dif.item.parms) <- all_items_names_dif

  #order base and dif item parms
  final_int_thr_base <- final$base.item.parms[grep(paste0(c(".int",".thr"),collapse = "|"),rownames(final$base.item.parms)),pen]
  final_slp_base <- final$base.item.parms[grep(".slp",rownames(final$base.item.parms)),pen]
  final_res_base <- final$base.item.parms[grep(".res",rownames(final$base.item.parms)),pen]
  final_names_base <- names(c(final_int_thr_base,final_slp_base,final_res_base))
  final$base.item.parms[,pen] <-  matrix(c(final_int_thr_base,final_slp_base,final_res_base),ncol=1)
  rownames(final$base.item.parms) <- final_names_base

  final_int_dif <- final$dif.item.parms[grep(".int",rownames(final$dif.item.parms)),pen]
  final_slp_dif <- final$dif.item.parms[grep(".slp",rownames(final$dif.item.parms)),pen]
  final_res_dif <- final$dif.item.parms[grep(".res",rownames(final$dif.item.parms)),pen]
  final_names_dif <- names(c(final_int_dif,final_slp_dif,final_res_dif))
  final$dif.item.parms[,pen] <-  matrix(c(final_int_dif,final_slp_dif,final_res_dif),ncol=1)
  rownames(final$dif.item.parms) <- final_names_dif

  #stop if there is a large change in DIF parameters
  if(pen > 1){
    second_last <- sum(final$dif.item.parms[-1,pen-1] == 0)
    last <- sum(final$dif.item.parms[-1,pen] == 0)
    if((second_last - last) > (num_predictors*num_items)){
      print(final)
      stop(paste0("Large increase in the number of DIF parameters from iteration ",pen-1," to ",pen,".\n  Two Options:\n  1. Provide smaller differences between lambda values.\n  2. Provide anchor item(s)."), call. = FALSE)
    }
  }

  #print information about optimization
  cat('\r',sprintf("Models Completed: %d of %d  Iteration: %d  Change: %d              ", pen, length(lambda), 0, 0))
  utils::flush.console()

  return(final)

}
