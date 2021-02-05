#' Pre-process data.
#'
#' @param item.data Matrix or data frame of item responses.
#' @param pred.data Matrix or data frame of DIF and/or impact predictors.
#' @param item.type Character value or vector indicating the item response
#' distributions.
#' @param pen.type Character indicating type of penalty.
#' @param tau Optional numeric vector of tau values.
#' @param num.tau Numeric indicating number of tau values to run Reg-DIF on.
#' @param anchor Optional numeric value or vector indicating which item
#' response(s) are anchors (e.g., \code{anchor = 1}).
#' @param stdz Logical value indicating whether to standardize DIF and
#' impact predictors for regularization.
#' @param control Optional list of additional model specification and
#' optimization parameters.
#' @param call Defined from regDIF.
#'
#' @keywords internal
#'
preprocess <-
  function(item.data,
           pred.data,
           item.type,
           pen.type,
           tau,
           num.tau,
           anchor,
           stdz,
           control,
           call){

  # Control parameters.
  final_control <- list(impact.mean.data = pred.data,
                        impact.var.data = pred.data,
                        tol = 10^-5,
                        maxit = 5000,
                        adapt.quad = FALSE,
                        num.quad = 21,
                        optim.method = "multi",
                        start.values = list())
  if(length(control) > 0) final_control[names(control)] <- control

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
  if(final_control$adapt.quad == T) {
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

  # Speed up computation.
  item_data <- as.matrix(sapply(item.data,as.numeric))
  pred_data <- as.matrix(sapply(pred.data,as.numeric))
  mean_predictors <-
    as.matrix(sapply(final_control$impact.mean.data,as.numeric))
  var_predictors <-
    as.matrix(sapply(final_control$impact.var.data,as.numeric))
  samp_size <- dim(item_data)[1]
  num_items <- dim(item_data)[2]
  num_predictors <- dim(pred_data)[2]

  # Get muliple characters of item.type for number of items.
  if(length(item.type) == 1 | is.null(item.type)) {
    if(is.null(item.type)) item.type <- "2PL"
    item_type <- rep(item.type, num_items)
  }

  # Get item response types.
  num_responses <- rep(1,num_items)
  cat_items <- item_type == "Rasch" |
    item_type == "2PL" |
    item_type == "Graded"
  if(any(cat_items)) {
    item_data[,which(cat_items)] <-
      apply(item.data[,which(cat_items)],
            2,
            function(x) as.numeric(as.factor(x)))
    num_responses[which(cat_items)] <-
      apply(item_data[,which(cat_items)],
            2,
            function(x) length(unique(na.omit(x))))
  }

  # Update number of quad pts for ordered categorical or guassian items.
  if(any(num_responses != 2) && final_control$num.quad == 21) {
    final_control$num.quad <- 51
  }

  # Define fixed quadrature points.
  theta <- seq(-6, 6, length.out = final_control$num.quad)

  if(final_control$adapt.quad == F && final_control$num.quad < 21) {
    warning(paste0("When using fixed quadrature, greater than 20 points for ",
                   "binary item responses or 50 points for ordered ",
                   "responses is recommended to yield precise estimates."))
  }

  # Standardize predictors.
  if(stdz == TRUE){
    pred_data <- scale(pred_data)
    mean_predictors <- scale(mean_predictors)
    var_predictors <- scale(var_predictors)
  }

  # Penalty type.
  if(is.null(pen.type)) {
    pen_type <- "lasso"
  } else {
    pen_type <- pen.type
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
      names(p[[item]]) <- c(paste0('c0_item',item,"_int",
                                   1:(num_responses[item]-1)),
                            paste0('a0_item',item,"_"),
                            paste0('c1_item',item,"_cov",1:num_predictors),
                            paste0('a1_item',item,"_cov",1:num_predictors))
    } else if(num_responses[item] == 2) {
      p[[item]] <- c(0,
                     1,
                     rep(0, num_predictors),
                     rep(0, num_predictors))
      names(p[[item]]) <- c(paste0('c0_item',item,"_int1"),
                            paste0('a0_item',item,"_"),
                            paste0('c1_item',item,"_cov",1:num_predictors),
                            paste0('a1_item',item,"_cov",1:num_predictors))
    } else if(num_responses[item] == 1) {
      p[[item]] <- c(mean(item_data[,item]),
                     sqrt(.5*var(item_data[,item])),
                     rep(0, num_predictors),
                     rep(0, num_predictors),
                     .5*var(item_data[,item]),
                     rep(0, num_predictors))
      names(p[[item]]) <- c(paste0('c0_item',item,"_int"),
                            paste0('a0_item',item, "_"),
                            paste0('c1_item',item,"_cov",1:num_predictors),
                            paste0('a1_item',item,"_cov",1:num_predictors),
                            paste0('s0_item',item, "_"),
                            paste0('s1_item',item,"_cov",1:num_predictors))
    }
  }

  p[[(num_items+1)]] <- rep(0,ncol(mean_predictors))
  p[[(num_items+2)]] <- rep(0,ncol(var_predictors))
  names(p[[(num_items+1)]]) <- paste0(rep(paste0('mean',
                                                 1:ncol(mean_predictors))))
  names(p[[(num_items+2)]]) <- paste0(rep(paste0('var',
                                                 1:ncol(var_predictors))))

  # Update starting values if provided by the user.
  if(length(final_control$start.values) > 0) {
    for(parm in 1:ncol(mean_predictors)) {
      p[[(num_items+1)]][[parm]] <- final_control$start.values$impact[parm]
    }
    for(parm in 1:ncol(var_predictors)) {
      p[[(num_items+2)]][[parm]] <-
        final_control$start.values$impact[ncol(mean_predictors)+parm]
    }
    for(item in 1:num_items) {
      p[[item]][[1]] <- final_control$start.values$base[item]
      p[[item]][[2]] <- final_control$start.values$base[num_items+item]
      for(cov in 1:num_predictors) {
        p[[item]][[2+cov]] <-
          final_control$start.values$dif[((item-1)*num_predictors)+cov]
        p[[item]][[length(p[[item]])-num_predictors+cov]] <-
          final_control$start.values$dif[(num_items*num_predictors)+((item-1)*num_predictors)+cov]
      }
    }
  }

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
  final_length <- ifelse(is.null(tau),num.tau,length(tau))
  final <- list(tau_vec = rep(NA,final_length),
                aic = rep(NA,final_length),
                bic = rep(NA,final_length),
                impact = matrix(NA,ncol=final_length,
                                         nrow=(ncol(mean_predictors)+
                                                 ncol(var_predictors))),
                base = matrix(NA,ncol=final_length,
                              nrow=num_base_parms),
                dif = matrix(NA,ncol=final_length,
                             nrow=num_dif_parms),
                em_history = lapply(1:num.tau,
                                    function(x) {
                                      mat <- matrix(0,
                                                    ncol=1,
                                                    nrow=length(unlist(p))+1)
                                      rownames(mat) <- c(names(unlist(p)),
                                                         "observed_ll")
                                      return(mat)
                                    }),
                log_like = NA,
                complete_ll_info <- list(),
                data = vector("list",2),
                call = call)
  return(list(p = p,
              final = final,
              item_data = item_data,
              pred_data = pred_data,
              mean_predictors = mean_predictors,
              var_predictors = var_predictors,
              item_type = item_type,
              theta = theta,
              final_control = final_control,
              num_responses = num_responses,
              num_predictors = num_predictors,
              samp_size = samp_size,
              num_items = num_items,
              pen_type = pen_type,
              tau_vec = tau_vec,
              num_tau = num.tau,
              id_tau = id_tau,
              num_quad = final_control$num.quad,
              adapt_quad = final_control$adapt.quad,
              optim_method = final_control$optim.method,
              em_history = final$em_history))

}
