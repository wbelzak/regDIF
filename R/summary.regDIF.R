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
      sum_results <- data.frame(
        "Lambda" = object$lambda[which.min(object$aic)],
        "AIC" = object$aic[which.min(object$aic)]
      )
      impact <- object$impact.lv.parms[,which.min(object$aic)]
      base <- object$base.item.parms[,which.min(object$aic)]
      dif <- object$dif.item.parms[,which.min(object$aic)]
    } else if(method == "bic"){
      sum_results <- data.frame(
        "Lambda" = object$lambda[which.min(object$bic)],
        "BIC" = object$bic[which.min(object$bic)]
      )
      impact <- object$impact.lv.parms[,which.min(object$bic)]
      base <- object$base.item.parms[,which.min(object$bic)]
      dif <- object$dif.item.parms[,which.min(object$bic)]
    }

    #get coef results for both min.lambda aic and bic

    impact <- rbind(impact[1:(length(impact)/2)],impact[-(1:(length(impact)/2))])
    colnames(impact) <- sub('.*mean.', '', colnames(impact))
    rownames(impact) <- c('Mean','Variance')

    base <- cbind(base[1:(length(base)/2)],base[-(1:(length(base)/2))])

    dif.int <- dif[1:(length(dif)/2)]
    dif.int <- split(dif.int,ceiling(seq_along(dif.int)/(nrow(base)/2)))
    dif.int <- lapply(dif.int,t)
    dif.int <- do.call(rbind,dif.int)
    dif.slp <- dif[-(1:(length(dif)/2))]
    dif.slp <- split(dif.slp,ceiling(seq_along(dif.slp)/(nrow(base)/2)))
    dif.slp <- lapply(dif.slp,t)
    dif.slp <- do.call(rbind,dif.slp)

    item.parms <- cbind(base[,1],dif.int,base[,2],dif.slp)
    colnames(item.parms) <- c("int.base",paste0("int.",colnames(impact)),"slp.base",paste0("slp.",colnames(impact)))
    rownames(item.parms) <- sub("\\..*", "", rownames(item.parms))

    rownames(sum_results) <- "Minimum"

    cat("\nregDIF Results:\n\n")
    ## print the results table
    print(sum_results)
    cat("\n")
    print(impact)
    cat("\n")
    print(item.parms)
  }
