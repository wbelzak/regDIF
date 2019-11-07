#' Regularized Differential Item Functioning
#'
#' Regularization of DIF parameters in item response theory (IRT) and moderated nonlinear factor analysis (MNLFA) models.
#'
#' @param data A matrix or data frame. Currently supports dichotomously (0-1) scored items only.
#' @param covariates A matrix or data frame of DIF covariates. Supports both categorical and continuous covariates.
#' @param tau A numeric vector of tuning parameters.
#' @param anchor A number or numeric vector indicating which items are anchors. Currently at least one item must be specified as an anchor, although more anchors may be specified. Default is item \code{1}.
#' @param quadpts The number of quadrature points to approximate the latent variable. More points lead to more precise estimates but slower run time. Default is \code{15} quadrature points.
#' @param standardize A number indicating which covariate to standardize. It is recommended that all covariates are on the same scale or standardized (i.e., Normal(0,1)). Default is \code{0} (no covariates are standardized).
#'
#' @return Function returns an object of class \code{SingleGroupClass}
#'
#' @examples
#' \dontrun{
#'
#' library(regDIF)
#' data <- ida[,1:6]
#' dif.covariates <- ida[,7:9]
#' tau <- seq(1:100)
#' model <- regDIF(data, dif.covariates, tau)
#' model
#'
#' }
#'
#' @import stats
#'
#' @export

regDIF <- function(data, covariates, tau, anchor = 1, quadpts = 15, standardize = 0){

  #stop to prevent improper data
  if(any(!(data == 0 | data == 1), na.rm = TRUE)) stop("Some data are not dichotomously scored as 0 or 1.")

  #get latent variable values (i.e., predictor values) for quadrature and tracelines
  theta <- seq(-4, 4, length.out = quadpts)

  #speeds up computation
  data <- as.matrix(data)
  covariates <- as.matrix(covariates)
  num_items <- dim(data)[2]
  samp_size <- dim(data)[1]
  num_quadpts <- length(theta)
  num_covariates <- dim(covariates)[2]
  num_parms <- num_items*num_covariates

  #standardize data
  if(standardize > 0 & length(standardize) == 1){
    covariates[,standardize] <- scale(covariates[,standardize])
  } else if(length(standardize) > 1){
    covariates[,standardize] <- sapply(covariates[,standardize], function(x) scale(x))
  }


#intial parameter values: baseline intercepts, baseline slopes, dif intercepts, dif slopes, impact
p <- c(rep(0, num_items), rep(1, num_items), rep(0, num_items*num_covariates*2), rep(0, num_covariates*2))
names(p) <- c(paste0('c0_itm',1:num_items),
              paste0('a0_itm',1:num_items),
              paste0(rep(paste0('c1_itm',1:num_items,'_cov'), each = num_covariates),rep(1:num_covariates, num_items)),
              paste0(rep(paste0('a1_itm',1:num_items,'_cov'), each = num_covariates),rep(1:num_covariates, num_items)),
              paste0(rep(paste0('g0_cov',1:num_covariates))),
              paste0(rep(paste0('b0_cov',1:num_covariates))))
final <- rep(list(rep(list(NA),2)),length(tau))

#Start Reg-DIF, loop through tuning parameters
for(t in 1:length(tau)){

  #Maximization/minimization routine: initialize the item parameters (starting values)
  if(t > 1){
    p_dif <- p[-grep(c('0'),names(p))]
    p_dif[p_dif < .001 & p_dif > -.001] <- 0
    p <- replace(p,names(p_dif),p_dif)
  }

  #lastp is the previous parameter estimates, in addition to maximization/minimization settings
  lastp <- p
  eps <- Inf
  tol = 10^-4
  maxit = 500
  iter = 0

  #loop until convergence or maximum number of iterations
  while(eps > tol & iter < maxit){

    #E-step:
    #evaluate E-step with current parameter estimates p
    rlist <- Estep.2pl(p,data,covariates,theta)

    #M-step 1: (impact)
    p_impact <- c(p[grep("g0",names(p))],p[grep("b0",names(p))])
    nr <- rlist[[3]]
    fit_impact = optim(par=p_impact,fn=ll.2pl.impact,nr=nr,theta=theta,covariates=covariates,method="BFGS",control=list(maxit = maxit))

    p <- replace(p,names(fit_impact$par),fit_impact$par)

    #M-step 2: (DIF)
    #for each item (loop); maximizing i independent logistic regression log-likelihoods (Q functions), with the quadrature points serving as the predictor values
    for (item in 1:num_items) {

      r1 <- rlist[[1]][[item]]
      r0 <- rlist[[2]][[item]]

      p_fit <- c(p[paste0("c0_itm",item)],p[grep(paste0("c1_itm",item),names(p))],p[paste0("a0_itm",item)],p[grep(paste0("a1_itm",item),names(p))])

      if(any(item == anchor)){
        p_fit <- p_fit[-c(grep(paste0("c1_itm",anchor[grep(item,anchor)]),names(p_fit)), grep(paste0("a1_itm",anchor[grep(item,anchor)]),names(p_fit)))]
      }

      if(t == 1){
        p_active <- p_fit
      } else if(t > 1){
        p_active <- p_fit[p_fit != 0]
      }

      fit_dif = optim(par=p_active,fn=ll.2pl.dif,r1=r1,r0=r0,theta=theta,covariates=covariates,tau=tau[t],method="BFGS",control=list(maxit = maxit))

      p <- replace(p,names(fit_dif$par),fit_dif$par)

    }

    #Update and check for convergence:
    #calculate the difference in parameter estimates
    eps = sqrt(sum((p - lastp)^2))

    #update parameter list
    lastp <- p

    #update the iteration number
    iter = iter + 1
    if(iter == maxit) warning("Iteration limit reached without convergence")
    cat(sprintf("Iteration: %d  Change: %f\n", iter, eps)) #print information about optimization


  } #end EM once converged or reached iteration limit

  #organize parameters into presentable form
  parms_impact <- t(matrix(c(p[grep("g0",names(p))],p[grep("b0",names(p))])))
  colnames(parms_impact) <- c(paste0('g0_cov',1:num_covariates),paste0('b0_cov',1:num_covariates))
  rownames(parms_impact) <- 'impact'

  parms_dif <- cbind(p[grep("c0_itm",names(p))],
                     do.call(rbind, split(p[grep("c1_itm",names(p))],rep(1:num_items,each=num_covariates))),
                     p[grep("a0_itm",names(p))],
                     do.call(rbind, split(p[grep("a1_itm",names(p))],rep(1:num_items,each=num_covariates))))
  colnames(parms_dif) <- c('c0',paste0('c',1:num_covariates),'a0',paste0('a',1:num_covariates))
  rownames(parms_dif) <- colnames(data)


  final[[t]][[1]] <- round(parms_impact,3)
  final[[t]][[2]] <- round(parms_dif,3)

  if(t > 1){
  if(sum(final[[t]][[2]][,-c(grep("c0",colnames(final[[t]][[2]])),grep("a0",colnames(final[[t]][[2]])))], na.rm = TRUE) == 0) break
  }

} #end Reg-DIF


return(final)
}

