#' Coefficient function for regDIF function
#'
#' @param object Fitted regDIF model object.
#' @param lambda Optional character or numeric indicating the lambda(s) at which the model coefficients are returned. For character value, may be \code{"lambda.min"}, which returns model coefficients for the value of lambda at which the minimum fit statistic is identified. For numeric, the value(s) provided corresponds to the value(s) of lambda.
#' @param method Character value indicating the model fit statistic to be used for determining \code{"lambda.min"}. Default is \code{"bic"}. May also be \code{"aic"}.
#' @param ... Additional arguments to be passed through to \code{coef}.
#' @rdname coef.regDIF
#'
#' @export

coef.regDIF <-
  function(object, lambda = NULL, method = "bic", ...) {

      #which lambda/model to show
      if(is.null(lambda)){
        table <- list("lambda.parms" = object$lambda,
                      "impact.lv.parms" = object$impact.lv.parms,
                      "base.item.parms" = object$base.item.parms,
                      "dif.item.parms" = object$dif.item.parms)
      } else if(lambda == "lambda.min"){
        if(method == "aic"){
          table <- list("lambda.parms" = object$lambda[which.min(object$aic)],
                        "impact.lv.parms" = object$impact.lv.parms[,which.min(object$aic)],
                        "base.item.parms" = object$base.item.parms[,which.min(object$aic)],
                        "dif.item.parms" = object$dif.item.parms[,which.min(object$aic)])
        }else if(method == "bic"){
          table <- list("lambda.parms" = object$lambda[which.min(object$bic)],
                        "impact.lv.parms" = object$impact.lv.parms[,which.min(object$bic)],
                        "base.item.parms" = object$base.item.parms[,which.min(object$bic)],
                        "dif.item.parms" = object$dif.item.parms[,which.min(object$bic)])
        }
      } else if(is.numeric(lambda)) {
        table <- list("lambda.parms" = object$lambda[lambda],
                      "impact.lv.parms" = object$impact.lv.parms[,lambda],
                      "base.item.parms" = object$base.item.parms[,lambda],
                      "dif.item.parms" = object$dif.item.parms[,lambda])
      }


    return(table)
  }


