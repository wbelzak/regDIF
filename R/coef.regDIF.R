#' Coefficient function for regDIF function
#'
#' @param object Fitted regDIF model object.
#' @param lambda Optional character or numeric indicating the lambda(s) at which the model coefficients are returned. For character value, may be \code{"lambda.min"}, which returns model coefficients for the value of lambda at which the minimum fit statistic is identified. For numeric, the value(s) provided corresponds to the value(s) of lambda.
#' @param method Character value indicating the model fit statistic to be used for determining \code{"lambda.min"}. Default is \code{"bic"}. May also be \code{"aic"}.
#' @param ... Additional arguments to be passed through.
#' @rdname coef.regDIF
#'
#' @export

coef.regDIF <-
  function(object, lambda = NULL, method = "bic", ...) {
    #create table to display results
    if(is.null(lambda)){
      table <- list("Lambda" = object$Lambda,
                    "Impact" = object$Impact,
                    "DIF" = object$DIF)
    } else if(lambda == "lambda.min"){
      if(method == "aic"){
        table <- list("Lambda" = object$Lambda[which.min(object$AIC)],
                      "Impact" = object$Impact[which.min(object$AIC)],
                      "DIF" = object$DIF[which.min(object$AIC)])
      }else if(method == "bic"){
        table <- list("Lambda" = object$Lambda[which.min(object$BIC)],
                      "Impact" = object$Impact[which.min(object$BIC)],
                      "DIF" = object$DIF[which.min(object$BIC)])
      }
    } else if(is.numeric(lambda)) {
      table <- list("Lambda" = object$Lambda[lambda],
                    "Impact" = object$Impact[lambda],
                    "DIF" = object$DIF[lambda])
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
