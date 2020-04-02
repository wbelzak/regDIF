###############
# Postprocess #
###############

postprocess <-
  function(estimates,
           responses,
           predictors,
           theta,
           lambda,
           pen,
           anchor,
           final.control,
           final,
           samp_size,
           num_responses,
           num_items) {

  #get estimates
  elist <- estimates[[1]]
  p <- estimates[[2]]

  #get information criteria
  infocrit <- information_criteria(elist,p,theta,predictors,lambda,pen,samp_size,num_responses,num_items,final.control$num_quadpts)

  #Organize impact parameters
  parms_impact <- rbind(p[[num_items+1]],p[[num_items+2]])
  if(is.null(colnames(x)) | length(colnames(x)) == 0){
    impact_colnames <- c(paste0('cov',1:num_predictors))
  } else{
    impact_colnames <- colnames(x)
  }
  colnames(parms_impact) <- impact_colnames
  rownames(parms_impact) <- c('Mean','Variance')

  #Organize dif parameters
  p2 <- unlist(p)
  parms_dif <- cbind(p2[grep("c0_itm",names(p2))],
                     do.call(rbind, split(p2[grep("c1_itm",names(p2),fixed=T)],rep(1:num_items,each=num_predictors))),
                     p2[grep("a0_itm",names(p2))],
                     do.call(rbind, split(p2[grep("a1_itm",names(p2),fixed=T)],rep(1:num_items,each=num_predictors))))
  if(is.null(colnames(x)) | length(colnames(x)) == 0){
    dif_colnames_int <- paste0('int_cov',1:num_predictors)
    dif_colnames_slp <- paste0('slp_cov',1:num_predictors)
  } else{
    dif_colnames_int <- paste0('int_',colnames(x))
    dif_colnames_slp <- paste0('slp_',colnames(x))
  }
  if(is.null(rownames(y)) | length(rownames(y)) == 0){
    dif_rownames <- paste0(1:num_items)
  } else{
    dif_rownames <- colnames(y)
  }
  colnames(parms_dif) <- c('int_base',
                           dif_colnames_int,
                           'slp_base',
                           dif_colnames_slp)
  rownames(parms_dif) <- dif_rownames

  #assign output to final list
  final$Lambda[pen] <- lambda[pen]
  final$AIC[pen] <- round(infocrit[1],2)
  final$BIC[pen] <- round(infocrit[2],2)
  final$Impact[[pen]] <- round(parms_impact,2)
  final$DIF[[pen]] <- round(parms_dif,2)


  #stop if lambda is too small on first run (this leads to different results b/c of identification constraints on DIF parameters)
  if(is.null(anchor) & pen == 1 & sum(abs(p2[grep(paste0("cov"),names(p2))])) > 0){
    print(coef(final))
    stop("First Lambda value is too small.\n  Two Options:\n  1. Increase first Lambda value large enough to ensure all DIF parameters are removed from the model.\n  2. Provide anchor item(s).", call. = TRUE)
  }
  #stop if there is a large change in DIF parameters
  if(pen > 1){
    second_last <- sum(final$DIF[[pen-1]][,-grep("0",colnames(final$DIF[[pen-1]]))] == 0)
    last <- sum(final$DIF[[pen]][,-grep("0",colnames(final$DIF[[pen-1]]))] == 0)
    if((second_last - last) > (num_predictors*num_items)){
      print(final)
      stop(paste0("Large increase in the number of DIF parameters from iteration ",pen-1," to ",pen,".\n  Two Options:\n  1. Provide smaller differences between lambda values.\n  2. Provide anchor item(s)."), call. = TRUE)
    }
  }

  return(final)

}
