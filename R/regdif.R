#' Regularized Differential Item Functioning
#'
#' Regularization of DIF parameters in item response theory (IRT) and moderated nonlinear factor analysis (MNLFA) models.
#'
#' @param data A matrix or data frame of item responses. Currently supports dichotomously (0-1) scored items only.
#' @param covariates A matrix or data frame of DIF covariates. Supports both categorical and continuous covariates.
#' @param tau A numeric vector of tuning parameter values.
#' @param anchor A number or numeric vector indicating which items are anchors. Currently at least one item must be specified as an anchor, although more anchors may be specified. Default is item \code{1}.
#' @param quadpts The number of quadrature points to approximate the latent variable. More points lead to more precise estimates but slower run time. Default is \code{15} quadrature points.
#' @param standardize Logical (\code{TRUE}/\code{FALSE}) indicating whether to normalize the covariates. Default is \code{TRUE}.
#'
#' @return Function returns an object of class \code{SingleGroupClass}
#'
#' @examples
#' \dontrun{
#'
#' library(regDIF)
#' data <- ida[,1:6]
#' dif.covariates <- ida[,7:9]
#' tau <- seq(2,0,-.05)
#' model <- regDIF(data, dif.covariates, tau)
#' model
#'
#' }
#'
#' @import stats
#'
#' @export

regDIF <- function(data, covariates, tau, anchor = 1, quadpts = 15, standardize = TRUE){


  ##############
  # Preprocess #
  ##############

  #stop to prevent improper data
  if(any(!(data == 0 | data == 1), na.rm = TRUE)) stop("Data must be scored 0 for no/incorrect and 1 for yes/correct.", call. = TRUE)

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
  final <- rep(list(rep(list(NA),3)),length(tau))

  ##############################
  # Reg-DIF - Loop through tau #
  ##############################

  for(t in 1:length(tau)){

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
      item <- Mstep.2pl.dif(lv[[1]],rlist,theta,covariates,tau,t,maxit,num_items,samp_size,num_quadpts,num_covariates,anchor)
      p <- item


      #Update and check for convergence: Calculate the difference in parameter estimates from current to previous
      eps = sqrt(sum((p - lastp)^2))

      #Update parameter list
      lastp <- p

      #update the iteration number
      iter = iter + 1
      if(iter == maxit) warning("EM iteration limit reached without convergence")
      cat(sprintf("****EM Iteration: %d  Change: %f\n", iter, eps)) #print information about optimization

    } #End EM once converged or reached iteration limit

    itemtrace <- replicate(n=num_items, matrix(0,nrow=samp_size,ncol=num_quadpts), simplify = F)
    ll_dif <- ll_dif_pen <- replicate(n=num_items, 0, simplify = F)
    for (item in 1:num_items) { #loop through items
      p_active <- c(p[paste0("c0_itm",item,"_")],p[grep(paste0("c1_itm",item,"_"),names(p),fixed=T)],p[paste0("a0_itm",item,"_")],p[grep(paste0("a1_itm",item,"_"),names(p),fixed=T)])
      itemtrace[[item]] <- trace.line.pts(p_active,theta,covariates) #computing probability of endorsement for theta value using current estimate of a and b
      ll_dif[[item]] <- (-1)*(sum(rlist[[1]][[item]]*(log(itemtrace[[item]])), na.rm = TRUE) + sum(rlist[[2]][[item]]*(log(1.0-itemtrace[[item]])), na.rm = TRUE))
      ll_dif_pen[[item]] <- ll_dif[[item]] - (tau[t]*sum(c(abs(p_active[grep("c1_itm",names(p_active))]), abs(p_active[grep("a1_itm",names(p_active))])), na.rm = TRUE)) #Q function we want to minimize
    }

    ll <- lv[[2]] + sum(unlist(ll_dif_pen))
    bic <- length(p[p!=0])*log(nrow(data)) + 2*ll

    ###############
    # Postprocess #
    ###############

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


    final[[t]][[1]] <- round(parms_impact,3)
    final[[t]][[2]] <- round(parms_dif,3)
    final[[t]][[3]] <- round(bic,3)

    #Stop Reg-DIF if all DIF parameters are equal to 0
    # if(t > 1){
    # if(sum(final[[t]][[2]][,-c(grep("c0",colnames(final[[t]][[2]])),grep("a0",colnames(final[[t]][[2]])))], na.rm = TRUE) == 0) {
    #   message(paste0("All DIF covariates have been removed from the model. Reg-DIF has terminated with current tau value = ",tau[t])) & break}
    # }

  } #Terminate Reg-DIF


#Obtain final results
return(final)

}

