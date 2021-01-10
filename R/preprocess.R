#' Pre-process data.
#'
#' @param item.data Matrix or data frame of item responses.
#' @param pred.data Matrix or data frame of DIF and/or impact predictors.
#' @param item.type Character value or vector indicating the item response
#' distributions.
#' @param num.tau Numeric value of how many to tau values to fit.
#' @param tau Optional numeric vector of tau values.
#' @param anchor Optional numeric value or vector indicating which item
#' response(s) are anchors (e.g., \code{anchor = 1}).
#' @param stdz Logical value indicating whether to standardize DIF and
#' impact predictors for regularization.
#' @param num.quad Numeric value indicating the number of quadrature points.
#' @param control Optional list of additional model specification and
#' optimization parameters.
#' @param pen.type Character indicating type of penalty.
#' @param adapt.quad Logical determining whether to use adaptive quadrature.
#' @param call Defined from regDIF.
#'
#' @keywords internal
#'
preprocess <-
  function(item.data,
           pred.data,
           item.type,
           num.tau,
           tau,
           anchor,
           stdz,
           num.quad,
           control,
           pen.type,
           adapt.quad,
           call){

  # Control parameters.
  final.control <- list(impact.data = list(mean = NULL, var = NULL),
                        tol = 10^-5,
                        maxit = 5000)
  if(length(control) > 0) final.control[names(control)] <- control

  # Impact data (if different).
  if(is.null(control$impact.data$mean)) {
    mean_predictors <- pred.data
  } else {
    mean_predictors <- control$impact.data$mean
  }

  if(is.null(control$impact.data$var)) {
    var_predictors <- pred.data
  } else {
    var_predictors <- control$impact.data$var
  }

  # Pre-process warnings.
  if(!(any(item.type == "Rasch") ||
       any(item.type == "2PL") ||
       any(item.type == "Graded") ||
       any(is.null(item.type)))) {
    stop(paste0("Item response types must either be Rasch, 2PL, or Graded."))
  }
  if(any(tau < 0)) {
    stop("Tau values must be non-negative.", call. = TRUE)
  }
  if(length(tau) > 1 && all(diff(tau) >= 0)) {
    stop("Tau values must be in descending order (e.g., tau = c(-2,-1,0)).")
  }
  if(is.null(anchor) && length(tau) == 1) {
    if(tau == 0) {
      stop("Anchor item must be specified with tau = 0.", call. = TRUE)
    }
  }
  if(!is.null(anchor) && !is.numeric(anchor)) {
    stop("Anchor items must be numeric (e.g., anchor = 1).", call. = TRUE)
  }
  if(adapt.quad == T) {
    warning(paste0("Adaptive quadrature is not fully supported. Fixed-point ",
                   "quadrature is recommended at this time."))
  }

  # Define number of tau values.
  if(is.null(tau)){
    tau_vec <- 1e20
    id_tau <- TRUE
  } else {
    num.tau <- length(tau)
    tau_vec <- tau
    id_tau <- FALSE
  }

  # Speeds up computation.
  item.data <- as.matrix(sapply(item.data,as.numeric))
  pred.data <- as.matrix(sapply(pred.data,as.numeric))
  mean_predictors <- as.matrix(sapply(mean_predictors,as.numeric))
  var_predictors <- as.matrix(sapply(var_predictors,as.numeric))
  samp_size <- dim(item.data)[1]
  num_items <- dim(item.data)[2]
  num_predictors <- dim(pred.data)[2]

  # Get muliple characters of item.type for number of items.
  if(length(item.type) == 1 | is.null(item.type)) {
    if(is.null(item.type)) item.type <- "2PL"
    item.type <- rep(item.type, num_items)
  }

  # Get item response types.
  num_responses <- rep(1,num_items)
  cat_items <- item.type == "Rasch" |
    item.type == "2PL" |
    item.type == "Graded"
  if(any(cat_items)) {
    item.data[,which(cat_items)] <-
      apply(item.data[,which(cat_items)],
            2,
            function(x) as.numeric(as.factor(x)))
    num_responses[which(cat_items)] <-
      apply(item.data[,which(cat_items)],
            2,
            function(x) length(unique(na.omit(x))))
  }


  # Define fixed quadrature points.
  if(is.null(num.quad) && all(num_responses == 2)) {
    num.quad <- 21
  } else if(is.null(num.quad) && any(num_responses > 2)) {
    num.quad <- 51
  }
  theta <- seq(-6, 6, length.out = num.quad)

  if(adapt.quad == F && num.quad < 21) {
    warning(paste0("When using fixed quadrature, greater than 20 points for ",
                   "binary item responses or 50 points for ordered ",
                   "responses is recommended to yield precise estimates."))
  }

  # Standardize predictors.
  if(stdz == TRUE){
    pred.data <- scale(pred.data)
    mean_predictors <- scale(mean_predictors)
    var_predictors <- scale(var_predictors)
  }

  # Penalty type.
  if(is.null(pen.type)) pen.type <- "lasso"

  # Starting values.
  p <- replicate(n=num_items+2,list(NA),simplify=F)
  for(item in 1:num_items){

    # Different item response types.
    if(num_responses[item] > 2) {
      p[[item]] <- c(0,
                     seq(.25,1,length.out = num_responses[item]-2),
                     1,
                     rep(0, num_predictors),
                     rep(0, num_predictors))
      names(p[[item]]) <- c(paste0('c0_itm',item,"_int"),
                            paste0('c0_itm',item,"_thr",
                                   1:(num_responses[item]-2),"_"),
                            paste0('a0_itm',item,"_"),
                            paste0('c1_itm',item,"_cov",1:num_predictors),
                            paste0('a1_itm',item,"_cov",1:num_predictors))
    } else if(num_responses[item] == 2) {
      p[[item]] <- c(0,
                     1,
                     rep(0, num_predictors),
                     rep(0, num_predictors))
      names(p[[item]]) <- c(paste0('c0_itm',item,"_int"),
                            paste0('a0_itm',item,"_"),
                            paste0('c1_itm',item,"_cov",1:num_predictors),
                            paste0('a1_itm',item,"_cov",1:num_predictors))
    } else if(num_responses[item] == 1) {
      p[[item]] <- c(mean(item.data[,item]),
                     sqrt(.5*var(item.data[,item])),
                     rep(0, num_predictors),
                     rep(0, num_predictors),
                     .5*var(item.data[,item]),
                     rep(0, num_predictors))
      names(p[[item]]) <- c(paste0('c0_itm',item,"_int"),
                            paste0('a0_itm',item,"_"),
                            paste0('c1_itm',item,"_cov",1:num_predictors),
                            paste0('a1_itm',item,"_cov",1:num_predictors),
                            paste0('s0_itm',item,"_"),
                            paste0('s1_itm',item,"_cov",1:num_predictors))
    }
  }

  p[[(num_items+1)]] <- rep(0,ncol(mean_predictors))
  p[[(num_items+2)]] <- rep(0,ncol(var_predictors))
  names(p[[(num_items+1)]]) <- paste0(rep(paste0('g',1:ncol(mean_predictors))))
  names(p[[(num_items+2)]]) <- paste0(rep(paste0('b',1:ncol(var_predictors))))

  if(any(item.type == "cfa")){
    num_base_parms <- length(c(unlist(p)[grep('c0',names(unlist(p)))],
                               unlist(p)[grep('a0',names(unlist(p)))],
                               unlist(p)[grep('s0',names(unlist(p)))]))
    num_dif_parms <- length(c(unlist(p)[grep('c1',names(unlist(p)))],
                              unlist(p)[grep('a1',names(unlist(p)))],
                              unlist(p)[grep('s1',names(unlist(p)))]))
  } else{
    num_base_parms <- length(c(unlist(p)[grep('c0',names(unlist(p)))],
                               unlist(p)[grep('a0',names(unlist(p)))]))
    num_dif_parms <- length(c(unlist(p)[grep('c1',names(unlist(p)))],
                              unlist(p)[grep('a1',names(unlist(p)))]))
  }
  final.length <- ifelse(is.null(tau),num.tau,length(tau))
  final <- list(tau_vec = rep(NA,final.length),
                aic = rep(NA,final.length),
                bic = rep(NA,final.length),
                impact.lv.parms = matrix(NA,ncol=final.length,
                                         nrow=(ncol(mean_predictors)+
                                                 ncol(var_predictors))),
                base.item.parms = matrix(NA,ncol=final.length,
                                         nrow=num_base_parms),
                dif.item.parms = matrix(NA,ncol=final.length,
                                        nrow=num_dif_parms),
                call = call)

  return(list(p = p,
              final = final,
              item.data = item.data,
              pred.data = pred.data,
              mean_predictors = mean_predictors,
              var_predictors = var_predictors,
              item.type = item.type,
              theta = theta,
              final.control = final.control,
              num_responses = num_responses,
              num_predictors = num_predictors,
              samp_size = samp_size,
              num_items = num_items,
              pen.type = pen.type,
              tau_vec = tau_vec,
              num.tau = num.tau,
              id_tau = id_tau,
              num.quad = num.quad))

}
