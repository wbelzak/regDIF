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

  #get estimates and information criteria
  p <- estimates[[1]]
  infocrit <- estimates[[2]]

  #Organize impact parameters
  lv_parms <- c(lambda[pen],p[[num_items+1]],p[[num_items+2]])
  if(is.null(colnames(x)) | length(colnames(x)) == 0){
    lv_names <- c('lambda',paste0('lv.mean.covariate',1:num_predictors),paste0('var.covariate',1:num_predictors))
  } else{
    lv_names <- c('lambda',paste0('lv.mean.',colnames(x)),paste0('lv.var.',colnames(x)))
  }

  #Organize item baseline parameters
  p2 <- unlist(p)
  item_parms_base <- c(lambda[pen],p2[grep("c0_itm",names(p2))],p2[grep("a0_itm",names(p2))])
  if(is.null(colnames(x)) | length(colnames(x)) == 0){
    item_names_base <- c('lambda',paste0('item',1:num_items,".c0"),paste0('item',1:num_items,".a0"))
  } else{
    item_names_base <- c('lambda',paste0('item',1:num_items,".c0.base"),paste0('item',1:num_items,".a0.base"))
  }

  #Organize item dif parameters
  item_parms_dif <- c(lambda[pen],p2[grep("c1_itm",names(p2),fixed=T)],p2[grep("a1_itm",names(p2),fixed=T)])
  if(is.null(colnames(x)) | length(colnames(x)) == 0){
    item_names_dif <- c('lambda',paste0(rep(paste0('item',1:num_items),each = 3),'.c',1:num_predictors),paste0(rep(paste0('item',1:num_items),each = 3),'.a',1:num_predictors))
  } else{
    item_names_dif <- c('lambda',paste0(rep(paste0('item',1:num_items),each = 3),'.c',1:num_predictors,'.',colnames(x)),paste0(rep(paste0('item',1:num_items),each = 3),'.a',1:num_predictors,'.',colnames(x)))
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

  #stop if lambda is too small on first run (this leads to different results b/c of identification constraints on DIF parameters)
  if(is.null(anchor) & pen == 1 & sum(abs(p2[grep(paste0("cov"),names(p2))])) > 0 & alpha == 1){
    print(coef(final))
    stop("First Lambda value is too small.\n  Two Options:\n  1. Increase first Lambda value large enough to ensure all DIF parameters are removed from the model.\n  2. Provide anchor item(s).", call. = TRUE)
  }

  #stop if there is a large change in DIF parameters
  if(pen > 1){
    second_last <- sum(final$DIF[[pen-1]][,-c(1,2+num_predictors)] == 0)
    last <- sum(final$DIF[[pen]][,-c(1,2+num_predictors)] == 0)
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
