#' Regularized Differential Item Functioning
#'
#' Regularization of DIF parameters in item response theory (IRT) and moderated nonlinear factor analysis (MNLFA) models.
#'
#' @param data Matrix or data frame of item responses. Currently supports dichotomously (0-1) scored items only.
#' @param covariates Matrix or data frame of DIF covariates. Supports both categorical and continuous covariates.
#' @param penalty Numeric vector of non-zero tuning parameter values. Must be in descending order, from largest to smallest values.
#' @param standardize Logical value indicating whether to normalize the covariates. Default is \code{TRUE}.
#' @param anchor Optional numeric vector indicating which items are anchors (e.g., \code{anchor = 1}). Default is \code{NULL}, which means that at least one DIF parameter per covariate across all items will be fixed to zero to identify the model.
#' @param quadpts Number of quadrature points to approximate the latent variable. More points lead to more precise estimates but slower run time. Default is \code{15} quadrature points.
#'
#' @return Function returns an object of class \code{SingleGroupClass}
#'
#' @examples
#' \dontrun{
#'
#' library(regDIF)
#' data <- ida[,1:6]
#' covariates <- ida[,7:9]
#' penalty <- seq(1,0,-.05)
#' model <- regDIF(data, covariates, penalty)
#' model
#'
#' }
#'
#' @import stats
#'
#' @export

regDIF <- function(data, covariates, penalty, standardize = TRUE, anchor = NULL, quadpts = 15){


  ##############
  # Preprocess #
  ##############

  #data warnings
  if(any(!(data == 0 | data == 1), na.rm = TRUE)) stop("Data must be scored 1 for yes/correct and 0 for no/incorrect.", call. = TRUE)
  #penalty warnings
  if(any(penalty < 0)) stop("Penalty values must be non-negative.", call. = TRUE)
  if(length(penalty) > 1 & all(diff(penalty) >= 0)) stop("Penalty values must be in descending order (e.g., penalty = c(-2,-1,0)).")
  #anchor warnings
  if(is.null(anchor) & length(penalty) == 1){if(penalty == 0) stop("Anchor item must be specified without penalty.", call. = TRUE)}
  if(!is.null(anchor) & !is.numeric(anchor)) stop("Anchor items must be numeric (e.g., anchor = 1).", call. = TRUE)
  #quadrature warnings
  if(quadpts < 1) stop("The number of quadrature points must be greater than 1.\n  15 or more points may be required to get accurate estimates.")
  if(!is.numeric(quadpts)) stop("Quadrature points must be numeric (e.g., quadpts = 15).")

  #get latent variable values (i.e., predictor values) for quadrature and tracelines
  theta <- seq(-4, 4, length.out = quadpts)

  #speeds up computation
  data <- as.matrix(data)
  covariates <- as.matrix(covariates)
  num_items <- dim(data)[2]
  samp_size <- dim(data)[1]
  num_quadpts <- length(theta)
  num_covariates <- dim(covariates)[2]

  #standardize data
  if(standardize == TRUE){
    covariates <- scale(covariates)
  }

  ###################
  # Starting Values #
  ###################
  p <- c(rep(0, num_items), rep(1, num_items), rep(0, num_items*num_covariates*2), rep(0, num_covariates*2))
  names(p) <- c(paste0('c0_itm',1:num_items,"_"),
                paste0('a0_itm',1:num_items,"_"),
                paste0(rep(paste0('c1_itm',1:num_items,'_cov'), each = num_covariates),rep(1:num_covariates, num_items)),
                paste0(rep(paste0('a1_itm',1:num_items,'_cov'), each = num_covariates),rep(1:num_covariates, num_items)),
                paste0(rep(paste0('g0_cov',1:num_covariates))),
                paste0(rep(paste0('b0_cov',1:num_covariates))))
  final <- rep(list(list(impact = NA,dif = NA,bic = NA,penalty = NA)),length(penalty))

  ##################################
  # Reg-DIF - Loop through penalty #
  ##################################
  message(paste0('Running ',length(penalty),if(length(penalty) == 1){' model.'}else{' models.'}))
  pb = txtProgressBar(min = 0, max = length(penalty), initial = 0)

  for(t in 1:length(penalty)){

    #Maximization settings
    lastp <- p
    eps <- Inf
    tol = 10^-5
    maxit = 1000
    iter = 0

    ################
    # EM Algorithm #
    ################

    #loop until convergence or maximum number of iterations
    while(eps > tol & iter < maxit){

      #E-step: Evaluate Q function with current parameter estimates p
      rlist <- Estep.2pl(p,data,covariates,theta,t,num_items,samp_size,num_quadpts)

      #M-step 1: Optimize impact parameters
      lv <- Mstep.2pl.impact(p,rlist,theta,covariates,maxit,samp_size,num_quadpts)
      p <- replace(p,names(lv[[1]]),lv[[1]])

      #M-step 2: Optimize DIF parameters
      p <- Mstep.2pl.dif(lv[[1]],data,rlist,theta,covariates,penalty,t,maxit,num_items,samp_size,num_quadpts,num_covariates,anchor)

      #Update and check for convergence: Calculate the difference in parameter estimates from current to previous
      eps = sqrt(sum((p - lastp)^2))

      #Update parameter list
      lastp <- p

      #update the iteration number
      iter = iter + 1
      if(iter == maxit) warning("EM iteration limit reached without convergence")
      # cat(sprintf("EM Iteration: %d  Change: %f\n", iter, eps)) #print information about optimization

    } #End EM once converged or reached iteration limit

    ###############
    # Postprocess #
    ###############

    #obtain likelihood values for all items
    itemtrace <- replicate(n=num_items, matrix(0,nrow=samp_size,ncol=num_quadpts), simplify = F)
    ll_dif <- ll_dif_pen <- pen <- replicate(n=num_items, 0, simplify = F)
    for (item in 1:num_items) { #loop through items
      p_active <- c(p[paste0("c0_itm",item,"_")],p[grep(paste0("c1_itm",item,"_"),names(p),fixed=T)],p[paste0("a0_itm",item,"_")],p[grep(paste0("a1_itm",item,"_"),names(p),fixed=T)])
      itemtrace[[item]] <- trace.line.pts(p_active,theta,covariates) #computing probability of endorsement for theta value using current estimate of a and b
      ll_dif[[item]] <- (-1)*(sum(rlist[[1]][[item]]*(log(itemtrace[[item]])), na.rm = TRUE) + sum(rlist[[2]][[item]]*(log(1.0-itemtrace[[item]])), na.rm = TRUE))
      pen <- (penalty[t]*sum(c(abs(p_active[grep("c1_itm",names(p_active))]), abs(p_active[grep("a1_itm",names(p_active))])), na.rm = TRUE))
      ll_dif_pen[[item]] <- ll_dif[[item]] - pen #Q function we want to minimize
    }

    ll <- lv[[2]] + sum(unlist(ll_dif_pen))
    bic <- length(p[p!=0])*log(nrow(data)) + 2*ll

    #organize parameters into presentable form
    parms_impact <- t(matrix(c(p[grep("g0",names(p),fixed=T)],p[grep("b0",names(p),fixed=T)])))
    colnames(parms_impact) <- c(paste0('g0_cov',1:num_covariates),paste0('b0_cov',1:num_covariates))
    rownames(parms_impact) <- 'impact'

    parms_dif <- cbind(p[grep("c0_itm",names(p))],
                       do.call(rbind, split(p[grep("c1_itm",names(p),fixed=T)],rep(1:num_items,each=num_covariates))),
                       p[grep("a0_itm",names(p))],
                       do.call(rbind, split(p[grep("a1_itm",names(p),fixed=T)],rep(1:num_items,each=num_covariates))))
    colnames(parms_dif) <- c('c0',paste0('c',1:num_covariates),'a0',paste0('a',1:num_covariates))
    rownames(parms_dif) <- colnames(data)

    final[[t]][[1]] <- round(parms_impact,2)
    final[[t]][[2]] <- round(parms_dif,2)
    final[[t]][[3]] <- round(bic,2)
    final[[t]][[4]] <- penalty[t]

    #stop if penalty is too small on first run (this leads to different results b/c of identification constraints on DIF parameters)
    if(is.null(anchor) & t == 1 & sum(abs(p[-grep(paste0("0"),names(p))])) > 0){
      print(final[[t]])
      stop("First penalty value is too small.\n  Increase first penalty value large enough to ensure all DIF parameters are removed from the model.\n  Or provide anchor item(s).", call. = TRUE)
    }

    #stop if there is a large change in DIF parameters
    if(t > 1){
      second_last <- sum(final[[t-1]]$dif[,-grep("0",colnames(final[[t-1]]$dif))] == 0)
      last <- sum(final[[t]]$dif[,-grep("0",colnames(final[[t]]$dif))] == 0)
      if((second_last - last) > (num_covariates*num_items)){
        print(final)
        stop(paste0("Large increase in the number of DIF parameters from iteration ",t-1," to ",t,".\n  Provide smaller differences between penalty values.\n  Or provide anchor item(s)."), call. = TRUE)
      }
    }

  setTxtProgressBar(pb,t)
  } #Terminate Reg-DIF


#Obtain final results
message(paste0('\nFinished.'))
return(final)

}

