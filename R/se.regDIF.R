#' Standard Errors for regDIF Model(s)
#'
#' Obtain standard errors for regDIF model(s).
#'
#' @usage
#' se.regDIF(model.obj,
#'           se.type = "sem",
#'           item.data = NULL,
#'           pred.data = NULL,
#'           tau = NULL,
#'           ...)
#'
#' @param model.obj A regDIF fitted model object of class \code{regDIF}.
#' Upon designating "model.obj", the default is to obtain standard errors for
#' the best-fitting model according to the minimum BIC model.
#' @param se.type Character value indicating the method of computing standard
#' errors for a regDIF fitted model. Default is "sem", or the supplemental EM
#' algorithm. Other options are not currently supported.
#' @param item.data Optional data frame or matrix of item response data to fit
#' regDIF model(s) and obtain standard errors. Default is NULL when
#' \code{model.obj} is provided. If providing \code{item.data}, then
#' \code{pred.data} and \code{tau} must also be provided.
#' @param pred.data Optional data frame or matrix of predictor data to fit
#' regDIF model(s) and obtain standard errors. Default is NULL when
#' \code{model.obj} is provided. If providing \code{pred.data}, then
#' \code{item.data} and \code{tau} must also be provided.
#' @param tau Optional numeric or vector of tau values to fit regDIF model(s)
#' and obtain standard errors. If providing \code{tau}, then either a
#' \code{model.obj} with a corresponding \code{tau} value (i.e., already fitted
#' to \code{tau}) must be provided or \code{item.data} and \code{pred.data}
#' must be provided.
#' @param ... Additional arguments to pass to regDIF when \code{item.data},
#' \code{pred.data}, and \code{tau} are provided.
#'
#' @return Function returns an object of class \code{se.regDIF}
#'
#' @examples
#' \dontrun{
#'
#' library(regDIF)
#' head(ida)
#' item.data <- ida[,1:6]
#' pred.data <- ida[,7:9]
#' fit <- regDIF(item.data, pred.data, tau = .2)
#' se.fit <- se.regDIF(fit)
#' se.fit
#'
#' }
#'
#' @import stats utils
#' @useDynLib regDIF, .registration = TRUE
#'
#' @export
se.regDIF <- function(model.obj,
                      se.type = "sem",
                      item.data = NULL,
                      pred.data = NULL,
                      tau = NULL,
                      ...) {

  if(is.null(item.data) &&
     is.null(pred.data) &&
     is.null(tau)) {

    # Obtain EM history of minimum BIC model from regDIF model object.
    em_history <- model.obj$em_history[[which.min(model.obj$bic)]]

    # Obtain MLEs from EM history.
    mle <- em_history[,ncol(em_history)]

    # Cycle through parameter estimates.
    em_history_updated <- mle
    em_history_before <- em_history[,ncol(em_history)-10]
    for(parm in 1:length(mle)) {

      # Skip estimates that are 0.
      if(parm == 0) next

      em_history_updated[parm] <- em_history_before[parm]

    }

  }


}
