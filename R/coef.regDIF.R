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
#'
#' @rdname coef.regDIF
#'
#' @return \code{NULL}
#' @export

coef.regDIF <-
  function(object, tau = NULL, method = "bic", ...) {

      # Which tau/model to show.
      if(is.null(tau)){
        table <- list("tau" = object$tau_vec,
                      "impact" = object$impact,
                      "base" = object$base,
                      "dif" = object$dif)

      } else if(tau == "tau.min") {
        if(method == "aic") {
          table <- list("tau" =
                          object$tau_vec[which.min(object$aic)],
                        "impact" =
                          object$impact[,which.min(object$aic)],
                        "base" =
                          object$base[,which.min(object$aic)],
                        "dif" =
                          object$dif[,which.min(object$aic)])

        } else if(method == "bic") {
          table <- list("tau" =
                          object$tau_vec[which.min(object$bic)],
                        "impact" =
                          object$impact[,which.min(object$bic)],
                        "base" =
                          object$base[,which.min(object$bic)],
                        "dif" =
                          object$dif[,which.min(object$bic)])
        }

      } else if(is.numeric(tau)) {
        table <- list("tau" =
                        object$tau_vec[tau],
                      "impact" =
                        object$impact[,tau],
                      "base" =
                        object$base[,tau],
                      "dif" =
                        object$dif[,tau])
      }


    return(table)
  }


