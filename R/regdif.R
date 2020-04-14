#' Regularized Differential Item Functioning
#'
#' Performs regularization of DIF effects in item response theory and confirmatory factor analysis models via penalized expectation-maximization.
#'
#' @usage
#' regDIF(x,
#'        y,
#'        family = c("bernoulli","categorical","gaussian"),
#'        penalty = c("lasso","mcp"),
#'        nlambda = 100,
#'        lambda.max = 2,
#'        alpha = 1,
#'        gamma = 3,
#'        lambda = NULL,
#'        anchor = NULL,
#'        rasch = FALSE,
#'        standardize = TRUE,
#'        quadpts = 81,
#'        control = list())
#'
#' @param x Matrix or dataframe of DIF predictors.
#' @param y Matrix or dataframe of item responses. See below for supported distributions.
#' @param family Character value or vector indicating the item response distributions via \code{y}. For scales where item responses are of one type only, the user may input one character value indicating the type (e.g., \code{"categorical"}). For mixed item types, the user must specify a vector of characters in the order that corresponds to the response matrix via \code{y}; e.g., \code{c(rep("categorical",2)}\code{, "bernoulli"}\code{, rep("gaussian",3))}. Supports:
#' \itemize{
#'    \item{\code{"bernoulli"} - Bernoulli item response via logistic link function (i.e., 1PL or 2PL model, see rasch option below for 1PL). Must be numeric/integer (2 unique values), factor (2 levels), or logical.}
#'    \item{\code{"categorical"} - Categorical item response via ordered logistic link function (i.e., Graded Response Model). Must be numeric/integer or factor.}
#'    \item{\code{"gaussian"} - Gaussian item response via identity link function (i.e., Confirmatory Factor Analysis). Must be numeric/integer.}}
#' @param penalty Character value indicating the penalty function to use. Supports:
#' \itemize{
#'    \item{\code{"lasso"} - The least absolute selection and shrinkage operator (LASSO) which controls DIF selection through \eqn{\lambda} (lambda).}
#'    \item{\code{"mcp"} - The minimax concave penalty (MCP) which controls DIF selection through \eqn{\lambda} (lambda) and estimator bias through \eqn{\gamma} (gamma).}}
#' @param nlambda Numeric value indicating how many lambda values to fit. Default is 100.
#' @param lambda.max Numberic value indicating the maximum lambda parameter to use for internal construction of lambda vector. Default is 3. Must be large enough to shrink all DIF effects to zero to begin with.
#' @param alpha Numeric value indicating the alpha parameter in the elastic net penalty function. Alpha controls the degree to which LASSO or ridge is used during regularization. Default is 1, which is equivalent to LASSO. For ridge, set alpha to 0. NOTE: If using MCP penalty, alpha may not be exactly 0.
#' @param gamma Numeric value indicating the gamma parameter in the MCP function. Gamma controls the degree of tapering of DIF effects as lambda decreases. Larger gamma leads to faster tapering (less bias but possibly more unstable optimization), whereas smaller gamma leads to slower tapering (more bias but more stable optimization). Default is 3. Must be greater than 1.
#' @param lambda Optional numeric vector of lambda values \eqn{\ge} 0. If lambda is supplied, this overrides the automatic construction of lambda values via \code{nlambda}. Must be non-negative and in descending order, from largest to smallest values (e.g., \code{seq(1,0,-.01)}.
#' @param anchor Optional numeric value or vector indicating which item response(s) are anchors (e.g., \code{anchor = 1}). Default is \code{NULL}, meaning at least one DIF effect per covariate will be fixed to zero as lambda approaches 0 (required to identify the model).
#' @param rasch Logical value indicating whether to constrain item slopes to 1 (i.e., equal slopes). If \code{TRUE}, no slope DIF will be evaluated. Default is \code{FALSE}.
#' @param standardize Logical value indicating whether to standardize DIF covariates for regularization. Default is \code{TRUE}, as it is recommended that all covariates be on the same scale.
#' @param quadpts Numeric value indicating the number of quadrature points to be used in approximating the latent variable distribution during estimation. Default is \code{81}. Not recommended to be below 10 for accurate estimation.
#' @param control Optional list of optimization parameters. May be:
#' \describe{
#'    \item{tol}{Convergence threshold of EM algorithm. Default is \code{10^-5}.}
#'    \item{maxiter}{Maximum number of EM iterations. Default is \code{10000}.}}
#'
#' @return Function returns an object of class \code{regDIF}
#'
#' @examples
#' \dontrun{
#'
#' library(regDIF)
#' head(ida)
#' x <- ida[,7:9]
#' y <- ida[,1:6]
#' fit <- regDIF(x, y, family = "bernoulli", penalty = "lasso")
#' fit
#'
#' }
#'
#' @import stats utils
#' @importFrom Rcpp sourceCpp
#' @useDynLib regDIF, .registration = TRUE
#'
#' @export

regDIF <- function(x,
                   y,
                   family = c("bernoulli","categorical","gaussian"),
                   penalty = c("lasso","mcp"),
                   nlambda = 100,
                   lambda.max = 2,
                   alpha = 1,
                   gamma = 3,
                   lambda = NULL,
                   anchor = NULL,
                   rasch = FALSE,
                   standardize = TRUE,
                   quadpts = 81,
                   control = list()){

  # data <- read.table("C:\\Users\\wbelz\\Dropbox\\Will\\Research\\Small Sample DIF Paper\\Small Samples Categorical Paper\\simdata\\ss500ni6dm2pd2im1\\ss500ni6dm2pd2im1_1.dat")
  # y <- data[,9:14]
  # y[which(y[,3]<3),3] <- 0
  # y[which(y[,3]>=3),3] <- 1
  # y[which(data[,10]==3 & data[,7] > 1.5),2] <- 4
  # x <- data[,3]
  # family <- c('categorical','categorical','bernoulli','categorical','categorical','categorical');penalty <- 'lasso';nlambda <- 100;lambda.max <- 2;alpha <- 1;pen <- 1;gamma <- 3;lambda <- .5;anchor <- 1;rasch <- F;standardize <- T;quadpts <- 81;control <- list()

  # family <- 'bernoulli';penalty <- 'lasso';nlambda <- 100;lambda.max <- 2;alpha <- 1;pen <- 1;gamma <- 3;lambda <- .5;anchor <- 1;rasch <- F;standardize <- T;quadpts <- 81;control <- list()
  # library(lavaan)
  # y <- HolzingerSwineford1939[,7:12]
  # x <- HolzingerSwineford1939[,2]
  # y[,1] <- ifelse(y[,1]>5,2,1)
  # y[,2] <- as.numeric(cut(y[,2],seq(min(y[,2]),max(y[,2]),by=1.75)))
  # # # write.table(cbind(1:nrow(y),y,x),"C:\\Users\\wbelz\\Dropbox\\Will\\Research\\Regularized IRT\\continuous_data.dat",col.names = F,row.names = F)
  # family = 'gaussian'
  # penalty = "lasso"
  # anchor = 1
  # standardize = F
  # rasch = F
  # nlambda <- 100;lambda.max <- 2;alpha <- 1;pen <- 1;gamma <- 3;lambda <- .5;anchor <- 1;rasch <- F;standardize <- T;quadpts <- 81;control <- list()

  # data_ds <- read.table("C:\\Users\\wbelz\\Dropbox\\Will\\Research\\Reg-DIF\\dissertation\\analysis\\scripts\\empirical\\ds.dat")
  # colnames(data_ds) <- c("id", "school", "male", "ageyrs", "live35", "maturity", "agecent", "agecent2", "maleage", "maleage2", "ds1", "ds2", "ds3", "ds4", "ds5", "ds6", "ds7", "ds8")
  # x <- data_ds[,c(3,7:10)]
  # y <- data_ds[,11:18]
  # family <- 'bernoulli';penalty <- 'lasso';nlambda <- 100;lambda.max <- 2;alpha <- 1;pen <- 1;gamma <- 3;lambda <- c(2,1);anchor <- 1;rasch <- F;standardize <- T;quadpts <- 81;control <- list()

  #obtain larger lambda if necessary
  need_larger_lambda <- TRUE #only true to start
  lambda_times <- 0
  while(need_larger_lambda | lambda_times < 6){


    #increase lambda.max if all penalized parameters have NOT been removed from model (unless specifying anchor item)
    if(need_larger_lambda == FALSE) break
    if(lambda_times > 0) {
      lambda.max <- lambda.max*1.5
      lambda[1] <- lambda[1]*1.5
    }
    #if too many lambda.max values have been tried, stop.
    if(lambda_times == 5){
      print(coef(data_scrub$final))
      stop("lambda.max is too small.\n  Three possible solutions:\n  1. Increase lambda.max large enough to ensure all DIF parameters are removed from the model.\n  2. Standardize predictors if not already standardized.\n  3. Provide anchor item(s).", call. = TRUE)
    }

    #preprocess data
    call <- match.call()
    data_scrub <- preprocess(x,y,family,penalty,nlambda,lambda.max,lambda,anchor,rasch,standardize,quadpts,control,call)


    #Run Reg-DIF by looping through lambda
    for(pen in 1:length(data_scrub$lambda)){

      #obtain regDIF estimates
      estimates <- em_estimation(data_scrub$p,data_scrub$responses,data_scrub$predictors,data_scrub$theta,data_scrub$itemtypes,penalty,data_scrub$lambda,alpha,gamma,pen,anchor,rasch,data_scrub$final.control,data_scrub$samp_size,data_scrub$num_items,data_scrub$num_responses,data_scrub$num_predictors,data_scrub$num_quadpts)

      #stop if lambda.max is too small on first run
      p2 <- unlist(estimates[[1]])
      dif_parms <- p2[grep(paste0("cov"),names(p2))]
      if(any(data_scrub$itemtypes == "gaussian")) dif_parms <- dif_parms[-grep("s1",names(dif_parms))]
      if(is.null(anchor) & #no anchor
         pen == 1 & #first penalty value
         sum(abs(dif_parms)) > 0 & #not all DIF effects are zero
         alpha == 1 #alpha is 1 for lasso
         ){
        message("\nWarning: lambda.max or user-defined lambda value is too small to penalize all parameters to zero without anchor item. Automatically trying larger lambda.max or lambda value.")
        lambda_times <- lambda_times + 1
        break
      } else{
        need_larger_lambda <- FALSE
      }

      #postprocess data
      data_final <- postprocess(estimates,data_scrub$responses,data_scrub$predictors,y,x,data_scrub$theta,data_scrub$lambda,alpha,pen,anchor,data_scrub$final.control,data_scrub$final,data_scrub$samp_size,data_scrub$num_responses,data_scrub$num_predictors,data_scrub$num_items,data_scrub$num_quadpts)

      #update parameter estimates for next lambda value
      data_scrub$p <- estimates[[1]]
      data_scrub$final <- data_final
    }

  }

  #Obtain final results
  class(data_final) <- "regDIF"
  cat("\n")
  return(data_final)

}

