#' Coefficient function for regDIF function
#'
#' @param regDIF_object Fitted regDIF model object.
#' @param lambda Optional character or numeric indicating the lambda(s) at which the model coefficients are returned. For character value, may be \code{"lambda.min"}, which returns model coefficients for the value of lambda at which the minimum fit statistic is identified. For numeric, the value(s) provided corresponds to the value(s) of lambda.
#' @param method Character value indicating the model fit statistic to be used for determining \code{"lambda.min"}. Default is \code{"bic"}. May also be \code{"aic"}.
#'
#' @return \code{NULL}
#' @export

coef.regDIF <-
  function(regDIF_object, lambda = NULL, method = "bic") {
    #create table to display results
    if(is.null(lambda)){
      table <- list("Lambda" = regDIF_object$Lambda,
                    "Impact" = regDIF_object$Impact,
                    "DIF" = regDIF_object$DIF)
    } else if(lambda == "lambda.min"){
      if(method == "aic"){
        table <- list("Lambda" = regDIF_object$Lambda[which.min(regDIF_object$AIC)],
                      "Impact" = regDIF_object$Impact[which.min(regDIF_object$AIC)],
                      "DIF" = regDIF_object$DIF[which.min(regDIF_object$AIC)])
      }else if(method == "bic"){
        table <- list("Lambda" = regDIF_object$Lambda[which.min(regDIF_object$BIC)],
                      "Impact" = regDIF_object$Impact[which.min(regDIF_object$BIC)],
                      "DIF" = regDIF_object$DIF[which.min(regDIF_object$BIC)])
      }
    } else if(is.numeric(lambda)) {
      table <- list("Lambda" = regDIF_object$Lambda[lambda],
                    "Impact" = regDIF_object$Impact[lambda],
                    "DIF" = regDIF_object$DIF[lambda])
    }
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
