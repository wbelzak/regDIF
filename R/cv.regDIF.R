#' Cross-Validation for regDIF
#'
#' Use k-fold cross validation to identify best-fitting model over a grid of
#' tau values (i.e., regularization tuning parameter).
#'
#' @usage
#' cv.regDIF(item.data,
#'        predictor.data,
#'        ...,
#'        nfolds = 10,
#'        fold,
#'        seed)
#'
#' @param item.data Matrix or dataframe of item responses, as in regDIF.
#' @param predictor.data Matrix or dataframe of DIF and/or impact predictors, as
#' in regDIF.
#' @param ... Additional arguments to regDIF.
#' @param nfolds Numeric value indicating the number of cross-validation folds.
#' The default is 10.
#' @param fold The fold that each observation belongs to. The default is to
#' assign the observations randomly.
#' @param seed Numeric value representing the seed for the random number
#' generator. Used to obtain reproducible results.
#' @return Function returns an object of class \code{cv.regDIF}
#'
#' @examples
#' \dontrun{
#'
#' library(regDIF)
#' head(ida)
#' item.data <- ida[,1:6]
#' predictor.data <- ida[,7:9]
#' cv.fit <- cv.regDIF(item.data, predictor.data)
#' summary(cv.fit)
#'
#' }
#'
#' @import stats utils
#'
#' @export
cv.regDIF <- function(item.data,
                   predictor.data,
                   ...,
                   nfolds = 10,
                   fold = NULL,
                   seed = 100) {

}
