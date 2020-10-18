#' Print function for regDIF function
#'
#' @param x Fitted regDIF model object.
#' @param ... Additional arguments to be passed through \code{print}.
#'
#' @rdname print.regDIF
#'
#' @return \code{NULL}
#' @export

print.regDIF <-
  function(x, ...) {
    ## print to screen with line break
    cat("Call:\n")
    ## print the model formula we fit
    print(x$call)
    ## create table to display results
    table <- data.frame(
      "Lambda" = x$lambda,
      "AIC" = x$aic,
      "BIC" = x$bic
    )
    cat("\nRegDIF Results:\n")
    ## print the results table
    print(table)
  }
