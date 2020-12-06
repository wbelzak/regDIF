#' Pre-process data.
#'
#' @param item.data Matrix or dataframe of item responses.
#' @param Matrix or dataframe of DIF and/or impact predictors.
#' @param item.type Character value or vector indicating the item response
#' distributions.
#' @param num.tau Numeric value of how many to tau values to fit.
#' @param max.tau Numberic value indicating the maximum tau parameter.
#' @param tau Optional numeric vector of tau values.
#' @param anchor Optional numeric value or vector indicating which item
#' response(s) are anchors (e.g., \code{anchor = 1}).
#' @param impact.data Optional list of matrices or data frames with predictors
#' for mean and variance impact.
#' @param standardize Logical value indicating whether to standardize DIF and
#' impact covariates for regularization.
#' @param quadpts Numeric value indicating the number of quadrature points.
#' @param control Optional list of optimization parameters.
#' @param call Defined from regDIF.
#'
#' @keywords internal
#'
preprocess <-
  function(item.data,
           predictor.data,
           item.type,
           num.tau,
           max.tau,
           tau,
           anchor,
           impact.data,
           standardize,
           quadpts,
           control,
           call){

  # Impact data (if different).
  if(is.null(impact.data$mean)) {
    mean_predictors <- predictor.data
  } else {
    mean_predictors <- impact.data$mean
  }

  if(is.null(impact.data$var)) {
    var_predictors <- predictor.data
  } else {
    var_predictors <- impact.data$var
  }

  # Pre-process warnings.
  # if(any(item.type == "cfa")) {
  #   stop(paste0("Continuous (Gaussian) item responses are not currently ",
  #               "supported."))
  # }
  if(!(any(item.type == "rasch") |
       any(item.type == "2pl") |
       any(item.type == "graded") |
       any(item.type == "cfa") |
       any(is.null(item.type)))) {
    stop(paste0("Item response types must either be rasch, 2pl, or ",
                "graded."))
  }
  if(any(tau < 0)) {
    stop("Tau values must be non-negative.", call. = TRUE)
  }
  if(length(tau) > 1 & all(diff(tau) >= 0)) {
    stop("Tau values must be in descending order (e.g., tau = c(-2,-1,0)).")
  }
  if(is.null(anchor) & length(tau) == 1) {
    if(tau == 0) {
      stop("Anchor item must be specified with tau = 0.", call. = TRUE)
    }

  }
  if(!is.null(anchor) & !is.numeric(anchor)) {
    stop("Anchor items must be numeric (e.g., anchor = 1).", call. = TRUE)
  }

  # Control parameters.
  final.control <- list(tol = 10^-5,
                        maxit = 10000)
  if(length(control) > 0) final.control[names(control)] <- control

  # Tau.
  if(is.null(tau)){
    tau <- seq(max.tau**(1/3),0,length.out = num.tau)**3
  }

  # Speeds up computation.
  item.data <- as.matrix(sapply(item.data,as.numeric))
  predictor.data <- as.matrix(sapply(predictor.data,as.numeric))
  mean_predictors <- as.matrix(sapply(mean_predictors,as.numeric))
  var_predictors <- as.matrix(sapply(var_predictors,as.numeric))
  samp_size <- dim(item.data)[1]
  num_items <- dim(item.data)[2]
  num_predictors <- dim(predictor.data)[2]

  # Get muliple characters of item.type for number of items.
  if(length(item.type) == 1 | is.null(item.type)) {
    if(is.null(item.type)) item.type <- "2pl"
    item.type <- rep(item.type, num_items)
  }

  # Get item response types.
  num_responses <- rep(1,num_items)
  cat_items <- item.type == "rasch" |
    item.type == "2pl" |
    item.type == "graded"
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

  # Standardize predictors.
  if(standardize == TRUE){
    predictor.data <- scale(predictor.data)
    mean_predictors <- scale(mean_predictors)
    var_predictors <- scale(var_predictors)
  }

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
  final <- list(tau = rep(NA,length(tau)),
                aic = rep(NA,length(tau)),
                bic = rep(NA,length(tau)),
                impact.lv.parms = matrix(NA,ncol=length(tau),
                                         nrow=(ncol(mean_predictors)+
                                                 ncol(var_predictors))),
                base.item.parms = matrix(NA,ncol=length(tau),
                                         nrow=num_base_parms),
                dif.item.parms = matrix(NA,ncol=length(tau),
                                        nrow=num_dif_parms),
                call = call)

  return(list(p = p,
              final = final,
              item.data = item.data,
              predictor.data = predictor.data,
              mean_predictors = mean_predictors,
              var_predictors = var_predictors,
              item.type = item.type,
              final.control = final.control,
              tau = tau,
              num_responses = num_responses,
              num_predictors = num_predictors,
              samp_size = samp_size,
              num_items = num_items))

}
