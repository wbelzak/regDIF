#' Summary function for regDIF function
#'
#' @param object Fitted regDIF model object.
#' @param method Fit statistic to use for displaying minimum tau model.
#' @param ... Additional arguments to be passed through \code{summary}.
#'
#' @rdname summary.regDIF
#'
#' @return \code{NULL}
#' @export

summary.regDIF <-
  function(object, method = "bic", ...) {
    # Print to screen with line break.
    cat("Call:\n")
    # Print the model formula we fit.
    print(object$call)
    # Create summary table to display results.
    if(method == "aic") {
      sum_results <- c(object$tau_vec[which.min(object$aic)],
                       object$aic[which.min(object$aic)])
      impact <- object$impact.lv.parms[,which.min(object$aic)]
      base <- object$base.item.parms[,which.min(object$aic)]
      dif <- object$dif.item.parms[,which.min(object$aic)]
      names(sum_results) <- c("tau","aic")

    } else if(method == "bic") {
      sum_results <- c(object$tau_vec[which.min(object$bic)],
                       object$bic[which.min(object$bic)])
      impact <- object$impact.lv.parms[,which.min(object$bic)]
      base <- object$base.item.parms[,which.min(object$bic)]
      dif <- object$dif.item.parms[,which.min(object$bic)]
      names(sum_results) <- c("tau","bic")

    }

    # Print the results table.
    cat(paste0("\nOptimal Model (out of ", length(object$tau_vec),"):\n"))
    print(sum_results)
    cat("\nLatent Variable Impact Parameters:\n")
    print(impact)
    cat("\nBase Item Parameters:\n")
    print(base)
    cat("\nDIF Item Parameters:\n")
    print(dif)
  }
