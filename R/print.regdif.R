#' Print function for regDIF function
#'
#' @param x regDIF model object to obtain basic results
#' @param ... Additional arguments to be passed through
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
      "Penalty" = x$Penalty,
      "BIC" = x$BIC
    )
    cat("\nRegDIF Results:\n")
    ## print the results table
    print(table)
  }
