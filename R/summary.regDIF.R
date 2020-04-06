#' Summary function for regDIF function
#'
#' @param object Fitted regDIF model object.
#' @param ... Additional arguments to be passed through \code{summary}.
#'
#' @rdname summary.regDIF
#'
#' @return \code{NULL}
#' @export

summary.regDIF <-
  function(object, ...) {
    ## print to screen with line break
    cat("Call:\n")
    ## print the model formula we fit
    print(object$call)
    ## create table to display results
    table <- data.frame(
      "Lambda" = object$lambda,
      "AIC" = object$aic,
      "BIC" = object$bic
    )
    cat("\nRegDIF Results:\n")
    ## print the results table
    print(table)
  }
