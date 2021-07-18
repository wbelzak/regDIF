#' Pre-process data.
#'
#' @param item.data Matrix or data frame of item responses.
#' @param pred.data Matrix or data frame of DIF and/or impact predictors.
#' @param prox.data Vector of observed proxy scores.
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
#' @return a \code{"list"} of default controls for \code{"em_estimation"}
#'
#' @keywords internal
#'
preprocess <-
  function(item.data,
           pred.data,
           prox.data,
           item.type,
           pen.type,
           tau,
           num.tau,
           anchor,
           stdz,
           control,
           call){

  # Remove observations with any NA.
  if(is.null(prox.data)) {
    combined.data <- cbind(pred.data, item.data)
  } else {
    combined.data <- cbind(pred.data, item.data, prox.data)
  }

  NA_cases <- apply(combined.data, 1, function(x) any(is.na(x)))
  if(sum(NA_cases) == nrow(combined.data)) {
    stop(paste0("No observations remain after performing listwise deletion. Consider removing ",
                "variables with the greatest amount of missingness."),
         call. = FALSE)
  }

  item.data <- item.data[!NA_cases,]
  pred.data <-
    if(is.null(dim(pred.data))) {
      pred.data[!NA_cases]
    } else {
      pred.data[!NA_cases,]
    }

  prox.data <- if(!is.null(prox.data)) prox.data[!NA_cases]

  # Control parameters.
  final_control <- list(impact.mean.data = pred.data,
                        impact.var.data = pred.data,
                        tol = 10^-5,
                        maxit = 2000,
                        adapt.quad = FALSE,
                        num.quad = 21,
                        int.limits = c(-6,6),
                        optim.method = "MNR",
                        start.values = list(),
                        parallel = list(FALSE,NULL))
  if(length(control) > 0) final_control[names(control)] <- control

  # Pre-process warnings.
  if(final_control$optim.method == "CD" && final_control$parallel[[1]]) {
    stop(paste0("Parallel computing is not supported for coordinate descent. Use \"UNR\" or ",
                "\"MNR\" for binary item responses or \"UNR\" for categorical item responses."))
  }
  if(!(any(item.type == "rasch") ||
       any(item.type == "2pl") ||
       any(item.type == "graded") ||
       any(item.type == "cfa") ||
       any(is.null(item.type)))) {
    stop("Item response types must be \"rasch\", \"2pl\", \"graded\", or \"cfa\".",
         call. = FALSE)
  }
  if(any(tau < 0)) {
    stop("Tau values must be non-negative.", call. = FALSE)
  }
  if(length(tau) > 1 && all(diff(tau) >= 0)) {
    stop("Tau values must be in descending order (e.g., tau = c(-2,-1,0)).", call. = FALSE)
  }
  if(is.null(anchor) && length(tau) == 1) {
    if(tau == 0) {
      stop("Anchor item must be specified with tau = 0.", call. = FALSE)
    }
  }
  if(!is.null(anchor) && !is.numeric(anchor)) {
    stop("Anchor items must be numeric (e.g., anchor = 1).", call. = FALSE)
  }
  if(!is.null(prox.data) && final_control$optim.method == "UNR") {
    stop("Coordinate descent is not yet supported when using observed proxy scores.", call. = FALSE)
  }
  if(final_control$adapt.quad == T) {
    warning(paste0("Adaptive quadrature is not fully supported. Fixed-point ",
                   "quadrature is recommended at this time."), call. = FALSE, immediate. = TRUE)
  }
  if(any(NA_cases)) {
    warning(paste0("Removed observations with missing values (NA)."), call. = FALSE,
    immediate. = TRUE)
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
  prox_data <- if(!is.null(prox.data)) scale(as.matrix(as.numeric(prox.data)))
  mean_predictors <-
    as.matrix(sapply(final_control$impact.mean.data,as.numeric))
  var_predictors <-
    as.matrix(sapply(final_control$impact.var.data,as.numeric))

  # Remove any variables with no variance.
  item_data_no_var <- apply(item_data, 2, var) == 0
  pred_data_no_var <- apply(pred_data, 2, var) == 0

  item_data <- item_data[, !item_data_no_var]
  pred_data <- pred_data[, !pred_data_no_var, drop = FALSE]

  if(any(item_data_no_var) || any(pred_data_no_var)) {
    warning(paste0("Removing the following variables from analysis because they have no variance ",
                   "(before or after list-wise deletion):"),
            paste0("\n"),
            paste0(names(item.data)[item_data_no_var],
                   names(pred.data)[pred_data_no_var], collapse = " "),
            call. = FALSE, immediate. = TRUE)
  }

  if(ncol(as.matrix(final_control$impact.mean.data)) != ncol(pred_data)) {
    mean_predictors <- mean_predictors[,!(names(pred.data) %in% names(pred.data)[pred_data_no_var])]
  }

  if(ncol(as.matrix(final_control$impact.var.data)) != ncol(pred_data)) {
    var_predictors <- var_predictors[,!(names(pred.data) %in% names(pred.data)[pred_data_no_var])]
  }

  # Get dimensions of data.
  samp_size <- dim(item_data)[1]
  num_items <- dim(item_data)[2]
  num_predictors <- dim(pred_data)[2]

  # Determine initial number of responses
  num_responses <-
    as.vector(apply(item_data, 2, function(x) length(unique(na.omit(x)))))


  # Get multiple characters of item.type for number of items.
  if(is.null(item.type)) {
    item_type <- sapply(1:num_items, function(item) {
      if(num_responses[item] == 2) {
        tmp = '2pl'
      } else if (num_responses[item] > 2 && num_responses[item] < 7) {
        tmp = 'graded'
      } else if (num_responses[item] > 6) {
        tmp = 'cfa'
      }
      return(tmp)
      })
  } else if(length(item.type) == 1) {
    item_type = rep(item.type, num_items)
  }

  # Get item response types.
  cat_items <- item_type == 'rasch' |
                     item_type == '2pl' |
                     item_type == 'graded'
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

  if(any(num_responses != 2)) {
    if(is.null(control$optim.method)) {
      warning(paste0("These item response types are not yet supported with ",
                     "Multivariate Newton-Raphson (MNR). Using Univariate Newton-Raphson ",
                     "(UNR) instead."), call. = FALSE, immediate. = TRUE)
    }
    final_control$optim.method <- "UNR"
  }

  # Update number of quad pts for ordered categorical or guassian items.
  if(any(num_responses != 2) && final_control$num.quad == 21) {
    if(any(num_responses > 6)) {
      final_control$num.quad <- 101
    } else {
      final_control$num.quad <- 51
    }
  }

  # Define fixed quadrature points.
  theta <- seq(final_control$int.limits[1],
               final_control$int.limits[2], length.out = final_control$num.quad)

  if(final_control$adapt.quad == F && final_control$num.quad < 21) {
    warning(paste0("When using fixed quadrature, greater than 20 points for ",
                   "binary item responses or 50 points for ordered ",
                   "responses is recommended to yield precise estimates."), call. = FALSE,
            immediate. = TRUE)
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
    if(item_type[item] == "cfa") {
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
    } else if(num_responses[item] == 2) {
      p[[item]] <- c(0,
                     1,
                     rep(0, num_predictors),
                     rep(0, num_predictors))
      names(p[[item]]) <- c(paste0('c0_item',item,"_int1"),
                            paste0('a0_item',item,"_"),
                            paste0('c1_item',item,"_cov",1:num_predictors),
                            paste0('a1_item',item,"_cov",1:num_predictors))
    } else {
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

  if(any(item_type == "cfa")){
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
                eap = list(scores = matrix(NA,ncol=final_length,
                                    nrow=length(NA_cases)),
                           sd = matrix(NA,ncol=final_length,
                                       nrow=length(NA_cases))),
                estimator_history = lapply(1:num.tau,
                                    function(x) {
                                      mat <- matrix(0,
                                                    ncol=1,
                                                    nrow=length(unlist(p))+1)
                                      rownames(mat) <- c(names(unlist(p)),
                                                         "observed_ll")
                                      return(mat)
                                    }),
                log_lik = NA,
                complete_ll_info = list(),
                data = vector("list",3),
                call = call)

  return(list(p = p,
              final = final,
              item_data = item_data,
              pred_data = pred_data,
              prox_data = prox_data,
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
              estimator_history = final$estimator_history,
              estimator_limit = F,
              exit_code = 0,
              NA_cases = NA_cases))

}
