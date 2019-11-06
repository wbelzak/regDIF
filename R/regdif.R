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

  #trace lines for item response functions
  trace.line.pts <-
    function(p_active,theta,covariates) {
      active_cov_c <- as.matrix(covariates[,as.numeric(gsub(".*([0-9]+)", "\\1", names(p_active[grep("c1",names(p_active))])))])
      active_cov_a <- as.matrix(covariates[,as.numeric(gsub(".*([0-9]+)", "\\1", names(p_active[grep("a1",names(p_active))])))])

      traceline <- 1/(1+exp(-((p_active[grep("c0",names(p_active))] + active_cov_c %*% p_active[grep("c1",names(p_active))]) + (p_active[grep("a0",names(p_active))] + active_cov_a %*% p_active[grep("a1",names(p_active))])*theta))) #logistic CDF

      }

  #log-likelihood function for impact
  ll.2pl.impact <-
    function(p_impact,nr,theta,covariates) {
      alpha <- covariates %*% p_impact[grep("g0",names(p_impact))]
      phi <- exp(covariates %*% p_impact[grep("b0",names(p_impact))])

      theta_scores <- matrix(0,nrow=samp_size,ncol=num_quadpts)
      for(case in 1:samp_size){
        theta_scores[case,] <- dnorm(theta, mean = alpha[case], sd = sqrt(phi[case]))
      }

      ll_impact <- (-1)*(sum(nr*log(theta_scores)))
    }

  #log-likelihood function for dif (conditional expected complete data log-likelihood)
  ll.2pl.dif <-
    function(p_active,r1,r0,theta,covariates,tau) { #p2 is current estimate of parameters, r1 is the conditional expected proportion of individuals that endorse item i, r0 is the conditional expected proportion of individuals that do not endorse item i

      itemtrace <- matrix(0,nrow=samp_size,ncol=num_quadpts)
      for(i in 1:num_quadpts){
        itemtrace[,i] <- trace.line.pts(p_active,theta[i],covariates)
      }

      ll_dif <- (-1)*(sum(r1*(log(itemtrace)), na.rm = TRUE) + sum(r0*(log(1.0-itemtrace)), na.rm = TRUE))
      pen <- (-1)*(tau*sum(c(abs(p_active[grep("c1_itm",names(p_active))]), abs(p_active[grep("a1_itm",names(p_active))])), na.rm = TRUE))
      ll_dif_pen <- ll_dif - pen #Q function we want to minimize
    }


  #Evaluating the Q function
  Estep.2pl <-
    function(p,data,covariates,theta) { #p is parameters, data is response data, theta is values of latent variable

      #remove non-active parameters from E-step
      if(t > 1){
        p <- p[p != 0]
      }

      #make space for the trace lines and the E-tables
      itemtrace <- replicate(n=num_items, matrix(0,nrow=samp_size,ncol=num_quadpts), simplify = F) #nrows = # of items, ncols = # of theta values
      r1 <- replicate(n=num_items, matrix(0,nrow=samp_size,ncol=num_quadpts), simplify = F)
      r0 <- replicate(n=num_items, matrix(0,nrow=samp_size,ncol=num_quadpts), simplify = F)


      #compute the trace lines
      for (item in 1:num_items) { #loop through items
        p_active <- c(p[paste0("c0_itm",item)],p[grep(paste0("c1_itm",item),names(p))],p[paste0("a0_itm",item)],p[grep(paste0("a1_itm",item),names(p))])
        for(q in 1:num_quadpts){
          itemtrace[[item]][,q] <- trace.line.pts(p_active,theta[q],covariates) #computing probability of endorsement for theta value using current estimate of a and b
        }
      }

      #impact
      alpha <- covariates %*% p[grep("g0",names(p))]
      phi <- exp(covariates %*% p[grep("b0",names(p))])

      #loop over responses to get r1 and r0 tables (conditional expected proportions of individuals that endorse/not endorse each item)
      for(case in 1:samp_size) { #looping over samples

        #getting qaudrature point ordinates (probablity weights)
        posterior <- dnorm(theta, mean = alpha[case], sd = sqrt(phi[case]))


        #within each response pattern, loop over items and compute posterior probability of response pattern
        for(item in 1:num_items) {
          x <- I(data[case,item]) #get one cell within data, **x will either be TRUE or FALSE** (endorse or not endorse)
          if (is.na(x))
            posterior <- posterior #if missing (NA), posterior probability remains the same as guassian points
          else if (x == 1)
            posterior <- posterior*itemtrace[[item]][case,] #if TRUE, prior probability weight times probability of endorsement
          else
            posterior <- posterior*(1-itemtrace[[item]][case,]) #if FALSE, prior probability weight times probability of no endorsement
        }


        #normalize posterior to area equals number of persons with this response pattern
        expd <- sum(posterior, na.rm = TRUE)
        posterior <- posterior/expd

        #for individual i, add posterior to the r1 and r0 tables depending on response
        for(item in 1:num_items) { #within a person, loop over items
          x <- I(data[case,item]) #get one cell within data, **x will either be TRUE or FALSE** (endorse or not endorse)
          if (is.na(x)) {
            r1[[item]][case,] <- r1[[item]][case,] #if missing (NA), conditional expected proportion of endorsing item j and not endorsing item j remains the same
            r0[[item]][case,] <- r0[[item]][case,]
          } else if (x == 1) {
            r1[[item]][case,] <- r1[[item]][case,] + posterior #if TRUE, add posterior to conditional expected proportion of endorsing item j
          } else
            r0[[item]][case,] <- r0[[item]][case,] + posterior #if FALSE, add posterior to conditoinal expected proportion of not endorsing item j
        }


      } #end loop over N
      nr <- Reduce('+',r0) + Reduce('+',r1)


      #list of r tables to be used in Q function (ll.2pl function above)
      rlist <- list(r1,r0,nr)


    } #end of E-step



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

