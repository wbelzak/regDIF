#' Coefficient function for regDIF function
#'
#' @param object regDIF model object to obtain coefficient values.
#' @param ... Additional arguments to be passed through
#'
#' @return \code{NULL}
#' @export

coef.regDIF <-
  function(object, ...) {
    ## create table to display results
    table <- list("Impact" = object$Impact,
                  "DIF" = object$DIF)
    ## print the results table
    return(table)
  }


# data <- data.frame(cbind(fit$Penalty,fit$BIC))
# colnames(data) <- c('penalty','bic')
# ggplot(data, aes(x=penalty,y=bic)) +
#   geom_line() +
#   scale_y_continuous(position = "right", breaks = seq(5460,5580,20)) +
#   scale_x_reverse(breaks = seq(.7,0,-.1), limits=c(.7, 0)) +
#   labs(title = "Lasso Path", x = expression(tau), y = "BIC") +
#   theme(axis.title.x = element_text(size=16,vjust = -.5),
#         axis.title.y.right = element_text(size=14,vjust = 2),
#         axis.text=element_text(size=12),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.background = element_blank(),
#         plot.title = element_text(hjust = 0.5),
#         axis.line = element_line(colour = "black"))
