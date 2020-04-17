#' Summary function for regDIF function
#'
#' @param object Fitted regDIF model object.
#' @param method Fit statistic to use for displaying minimum lambda model.
#' @param ... Additional arguments to be passed through \code{summary}.
#'
#' @rdname summary.regDIF
#'
#' @return \code{NULL}
#' @export

summary.regDIF <-
  function(object, method = "bic", ...) {
    ## print to screen with line break
    cat("Call:\n")
    ## print the model formula we fit
    print(object$call)
    ## create summary table to display results
    if(method == "aic"){
      sum_results <- c(object$lambda[which.min(object$aic)],object$aic[which.min(object$aic)])
      impact <- object$impact.lv.parms[,which.min(object$aic)]
      base <- object$base.item.parms[,which.min(object$aic)]
      dif <- object$dif.item.parms[,which.min(object$aic)]
      names(sum_results) <- c("lambda","aic")
    } else if(method == "bic"){
      sum_results <- c(object$lambda[which.min(object$bic)],object$bic[which.min(object$bic)])
      impact <- object$impact.lv.parms[,which.min(object$bic)]
      base <- object$base.item.parms[,which.min(object$bic)]
      dif <- object$dif.item.parms[,which.min(object$bic)]
      names(sum_results) <- c("lambda","bic")
    }

    #get coef results for both min.lambda aic and bic

    # #impact
    # impact <- rbind(impact[1:(length(impact)/2)],impact[-(1:(length(impact)/2))])
    # colnames(impact) <- sub('.*mean.', '', colnames(impact))
    # rownames(impact) <- c('Mean','Variance')
    #
    # #base parameters
    # if(length(grep(".thr",names(base))) > 0){ #if threshold values are included
    #   int_parms_base <- base[grep(".int.base",names(base))]
    #   thr_parms_base <- base[grep(".thr",names(base))]
    #   slp_parms_base <- base[grep(".slp.base",names(base))]
    #
    #   int_names <- sub("\\..*","",names(int_parms_base))
    #   thr_names <- sub("\\..*","",names(thr_parms_base))
    #   thr_table <- table(match(thr_names,int_names))
    #   thr_table_names <- names(thr_table)
    #
    #   no_thr <- 1:length(int_names) %in% as.numeric(names(thr_table))
    #   thr_table <- c(thr_table,rep(NA,length(which(!no_thr))))
    #   names(thr_table) <- c(thr_table_names,which(!no_thr))
    #   thr_table <- thr_table[order(names(thr_table))]
    #
    #   parms_base <- cbind(int_parms_base,matrix(NA, nrow = length(int_parms_base), ncol = max(thr_table,na.rm=T)),slp_parms_base)
    #   for(j in 1:length(int_parms_base)){
    #     if(!is.na(thr_table[j])){
    #       k <- j+(thr_table[j]-1)
    #       parms_base[j,2:(1+thr_table[j])] <- thr_parms_base[j:k]
    #     } else{
    #       parms_base[j,2] <- NA
    #     }
    #   }
    #
    # }
    #
    # if(length(grep(".res.base",names(base))) > 0){ #if residual parameters are included
    #   int_parms_base <- base[grep(".int.base",names(base))]
    #   res_parms_base <- base[grep(".res.base",names(base))]
    #
    #   int_names <- sub("\\..*","",names(int_parms_base))
    #   res_names <- sub("\\..*","",names(res_parms_base))
    #   res_table <- table(match(res_names,int_names))
    #   res_table_names <- names(res_table)
    #
    #   no_res <- 1:length(int_names) %in% as.numeric(names(res_table))
    #   res_table <- c(res_table,rep(NA,length(which(!no_res))))
    #   names(res_table) <- c(res_table_names,which(!no_res))
    #   res_table <- res_table[order(names(res_table))]
    #
    #   parms_base <- cbind(parms_base,matrix(NA, nrow = length(int_parms_base), ncol = 1))
    #   for(j in 1:length(int_parms_base)){
    #     if(!is.na(thr_table[j])){
    #       k <- j+(thr_table[j]-1)
    #       parms_base[j,2:(1+thr_table[j])] <- thr_parms_base[j:k]
    #     } else{
    #       parms_base[j,2] <- NA
    #     }
    #   }
    # }
    #
    # dif.int <- dif[1:(length(dif)/2)]
    # dif.int <- split(dif.int,ceiling(seq_along(dif.int)/(length(dif)/2)))
    # dif.int <- lapply(dif.int,t)
    # dif.slp <- dif[-(1:(length(dif)/2))]
    # dif.slp <- split(dif.slp,ceiling(seq_along(dif.slp)/(length(dif)/2)))
    # dif.slp <- lapply(dif.slp,t)
    # if(length(dif)/2 == length(dif.int))
    # dif.int <- do.call(rbind,dif.int)
    # dif.slp <- do.call(rbind,dif.slp)
    #
    # item.parms <- cbind(parms_base[,1:(1+max(thr_table,na.rm=T))],dif.int,parms_base[,(2+max(thr_table,na.rm=T)):ncol(parms_base)],dif.slp)
    # colnames(item.parms) <- c("int.base",paste0("int.",colnames(impact)),"slp.base",paste0("slp.",colnames(impact)))
    # rownames(item.parms) <- sub("\\..*", "", rownames(item.parms))
    #

    ## print the results table
    cat(paste0("\nOptimal Model (out of ", length(object$lambda),"):\n"))
    print(sum_results)
    cat("\nLatent Variable Impact Parameters:\n")
    print(impact)
    cat("\nBase Item Parameters:\n")
    print(base)
    cat("\nDIF Item Parameters:\n")
    print(dif)
  }
