###############
# Postprocess #
###############

postprocess <-
  function(estimates,
           responses,
           predictors,
           y,
           x,
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

    # responses <- data_scrub$responses;predictors <- data_scrub$predictors;theta <- data_scrub$theta;itemtypes <- data_scrub$itemtypes;lambda <- data_scrub$lambda; final.control <- data_scrub$final.control;final <- data_scrub$final;samp_size <- data_scrub$samp_size;num_items <- data_scrub$num_items;num_responses <- data_scrub$num_responses;num_predictors <- data_scrub$num_predictors;num_quadpts <- data_scrub$num_quadpts

  #get estimates and information criteria
  p <- estimates[[1]]
  infocrit <- estimates[[2]]

  #Organize impact parameters
  lv_parms <- c(p[[num_items+1]],p[[num_items+2]])
  if(is.null(colnames(x)) | length(colnames(x)) == 0){
    lv_names <- c(paste0('lv.mean.covariate',1:num_predictors),paste0('var.covariate',1:num_predictors))
  } else{
    lv_names <- c(paste0('lv.mean.',colnames(x)),paste0('lv.var.',colnames(x)))
  }

  #Organize item baseline parameters
  p2 <- unlist(p)
  item_parms_base <- c(p2[grep("c0_itm",names(p2))],p2[grep("a0_itm",names(p2))])
  if(is.null(colnames(y)) | length(colnames(y)) == 0){
    item_names_base <- c(paste0('item',1:num_items,".int.base"),paste0('item',1:num_items,".slp.base"))
  } else{
    item_names_base <- c(paste0(colnames(y),".int.base"),paste0(colnames(y),".slp.base"))
  }

  #Organize item dif parameters
  item_parms_dif <- c(p2[grep("c1_itm",names(p2),fixed=T)],p2[grep("a1_itm",names(p2),fixed=T)])
  if(is.null(colnames(x)) | length(colnames(x)) == 0){
    item_names_dif <- c(paste0(rep(paste0('item',1:num_items),each = 3),'.int.cov',1:num_predictors),paste0(rep(paste0('item',1:num_items),each = 3),'.slp.cov',1:num_predictors))
  } else{
    item_names_dif <- c(paste0(rep(paste0('item',1:num_items),each = 3),'.int.',colnames(x)),paste0(rep(paste0('item',1:num_items),each = 3),'.slp.',colnames(x)))
  }



  #assign output to final list
  final$lambda[pen] <- lambda[pen]
  final$aic[pen] <- round(infocrit[1],3)
  final$bic[pen] <- round(infocrit[2],3)
  final$impact.lv.parms[,pen] <- round(lv_parms,3)
  final$base.item.parms[,pen] <- round(item_parms_base,3)
  final$dif.item.parms[,pen] <- round(item_parms_dif,3)
  rownames(final$impact.lv.parms) <- lv_names
  rownames(final$base.item.parms) <- item_names_base
  rownames(final$dif.item.parms) <- item_names_dif

  #stop if there is a large change in DIF parameters
  if(pen > 1){
    second_last <- sum(final$dif.item.parms[-1,pen-1] == 0)
    last <- sum(final$dif.item.parms[-1,pen] == 0)
    if((second_last - last) > (num_predictors*num_items)){
      print(final)
      stop(paste0("Large increase in the number of DIF parameters from iteration ",pen-1," to ",pen,".\n  Two Options:\n  1. Provide smaller differences between lambda values.\n  2. Provide anchor item(s)."), call. = TRUE)
    }
  }

  #print information about optimization
  cat('\r',sprintf("Models Completed: %d of %d  Iteration: %d  Change: %d              ", pen, length(lambda), 0, 0))
  utils::flush.console()

  return(final)

}
