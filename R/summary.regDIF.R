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
      dif <- object$dif[object$dif[,which.min(object$aic)] != 0, which.min(object$aic)]
      if(length(grep(".res.", names(dif))) != 0) dif <- dif[-grep(".res.", names(dif))]
      names(sum_results) <- c("tau","aic")

    } else if(method == "bic") {
      sum_results <- c(object$tau_vec[which.min(object$bic)],
                       object$bic[which.min(object$bic)])
      dif <- object$dif[object$dif[,which.min(object$bic)] != 0, which.min(object$bic)]
      if(length(grep(".res.", names(dif))) != 0) dif <- dif[-grep(".res.", names(dif))]
      names(sum_results) <- c("tau","bic")

    }

    # Print the results table.
    cat(paste0("\nOptimal model (out of ", length(object$tau_vec),"):\n"))
    print(sum_results)
    cat("\nNon-zero DIF effects:\n")
    if(length(dif) != 0) print(dif)
  }
