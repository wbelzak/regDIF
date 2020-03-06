#' Coefficient function for regDIF function
#'
#' @return \code{NULL}
#' @export

coef.regDIF <-
  function(x) {
    ## print to screen with line break
    cat("Call:\n")
    ## print the model formula we fit
    print(x$call)
    ## create table to display results
    table <- list("Impact" = x$Impact,
                  "DIF" = x$DIF)
    cat("\nRegDIF Results:\n")
    ## print the results table
    print(table)
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