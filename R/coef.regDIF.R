#' Coefficient function for regDIF function
#'
#' @param object Fitted regDIF model object.
#' @param tau Optional character or numeric indicating the tau(s) at
#' which the model coefficients are returned. For character value, may be
#' \code{"tau.min"}, which returns model coefficients for the value of tau
#' at which the minimum fit statistic is identified. For numeric, the value(s)
#' provided corresponds to the value(s) of tau.
#' @param method Character value indicating the model fit statistic to be used
#' for determining \code{"tau.min"}. Default is \code{"bic"}. May also be
#' \code{"aic"}.
#' @param ... Additional arguments to be passed through to \code{coef}.
#' @rdname coef.regDIF
#'
#' @export

coef.regDIF <-
  function(object, tau = NULL, method = "bic", ...) {

      # Which tau/model to show.
      if(is.null(tau)){
        table <- list("tau.parms" = object$tau,
                      "impact.lv.parms" = object$impact.lv.parms,
                      "base.item.parms" = object$base.item.parms,
                      "dif.item.parms" = object$dif.item.parms)

      } else if(tau == "tau.min") {
        if(method == "aic") {
          table <- list("tau.parms" =
                          object$tau[which.min(object$aic)],
                        "impact.lv.parms" =
                          object$impact.lv.parms[,which.min(object$aic)],
                        "base.item.parms" =
                          object$base.item.parms[,which.min(object$aic)],
                        "dif.item.parms" =
                          object$dif.item.parms[,which.min(object$aic)])

        } else if(method == "bic") {
          table <- list("tau.parms" =
                          object$tau[which.min(object$bic)],
                        "impact.lv.parms" =
                          object$impact.lv.parms[,which.min(object$bic)],
                        "base.item.parms" =
                          object$base.item.parms[,which.min(object$bic)],
                        "dif.item.parms" =
                          object$dif.item.parms[,which.min(object$bic)])
        }

      } else if(is.numeric(tau)) {
        table <- list("tau.parms" =
                        object$tau[tau],
                      "impact.lv.parms" =
                        object$impact.lv.parms[,tau],
                      "base.item.parms" =
                        object$base.item.parms[,tau],
                      "dif.item.parms" =
                        object$dif.item.parms[,tau])
      }


    return(table)
  }


