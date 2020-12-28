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
    # Print to screen with line break.
    cat("Call:\n")
    # Print the model formula we fit.
    print(x$call)
    # Create table to display results.
    table <- data.frame(
      "Tau" = x$tau_vec,
      "AIC" = x$aic,
      "BIC" = x$bic
    )
    cat("\nRegDIF Results:\n")
    # Print the results table.
    print(table)
  }
