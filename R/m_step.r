######################
# M-step 1: (Impact) #
######################

Mstep.2pl.impact <-
  function(p,rlist,theta,covariates,maxit,samp_size,num_quadpts){

    p_impact <- c(p[grep("g0",names(p))],p[grep("b0",names(p))])
    nr <- rlist[[3]]
    #M-step 1: (Impact)
    fit_impact = optim(par=p_impact,fn=ll.2pl.impact,nr=nr,theta=theta,covariates=covariates,samp_size=samp_size,num_quadpts=num_quadpts,method="BFGS",control=list(maxit = maxit))

    p <- replace(p,names(fit_impact$par),fit_impact$par)
    ll <- fit_impact$value
    return(list(p,ll))

  }


###################
# M-step 2: (DIF) #
###################
Mstep.2pl.dif <-
  function(p,rlist,theta,covariates,tau,t,maxit,num_items,samp_size,num_quadpts,anchor){
    # ll <- 0
    #for each item (loop); maximizing i independent logistic regression log-likelihoods (Q functions), with the quadrature points serving as the predictor values
    for (item in 1:num_items) {

      r1 <- rlist[[1]][[item]]
      r0 <- rlist[[2]][[item]]

      lastp <- p
      eps <- Inf
      tol = 10^-3
      maxit = 20
      iter = 0

      while(eps > tol & iter < maxit){

      anl_deriv <- d("c0",p,r1,r0,data,item.focus=item,theta,covariates,cov=NULL,num_items,samp_size,num_quadpts)
      p_active <- c(p[paste0("c0_itm",item,"_")],p[grep(paste0("c1_itm",item,"_"),names(p))],p[paste0("a0_itm",item,"_")],p[grep(paste0("a1_itm",item,"_"),names(p))])
      p_current <- p_active[paste0("c0_itm",item,"_")]
      p_new <- p_current - anl_deriv[[1]]/anl_deriv[[2]]
      p <- replace(p,names(p_new),p_new)

      anl_deriv <- d("a0",p,r1,r0,data,item.focus=item,theta,covariates,cov=NULL,num_items,samp_size,num_quadpts)
      p_active <- c(p[paste0("c0_itm",item,"_")],p[grep(paste0("c1_itm",item,"_"),names(p))],p[paste0("a0_itm",item,"_")],p[grep(paste0("a1_itm",item,"_"),names(p))])
      p_current <- p_active[paste0("a0_itm",item,"_")]
      p_new <- p_current - anl_deriv[[1]]/anl_deriv[[2]]
      p <- replace(p,names(p_new),p_new)

      if(any(item == anchor)){
        #Update and check for convergence: Calculate the difference in parameter estimates from current to previous
        eps = sqrt(sum((p - lastp)^2))

        #Update parameter list
        lastp <- p

        #update the iteration number
        iter = iter + 1
        if(iter == maxit) warning("Coordinate descent iteration limit reached without convergence")

        next
      }

      p_c1 <- p[grep(paste0("c1_itm",item),names(p))]
      p_a1 <- p[grep(paste0("a1_itm",item),names(p))]

      for(cov in 1:length(p_c1)){
        anl_deriv <- d("c1",p,r1,r0,data,item.focus=item,theta,covariates,cov,num_items,samp_size,num_quadpts)
        p_active <- c(p[paste0("c0_itm",item,"_")],p[grep(paste0("c1_itm",item,"_"),names(p))],p[paste0("a0_itm",item,"_")],p[grep(paste0("a1_itm",item,"_"),names(p))])
        p_current <- p[grep(paste0("c1_itm",item,"_cov",cov),names(p))]
        z <- (anl_deriv[[2]]*p_current - anl_deriv[[1]])/anl_deriv[[2]]
        p_new <- sign(z)*max(abs(z) - tau[t], 0)
        p <- replace(p,names(p_new),p_new)
      }

      for(cov in 1:length(p_a1)){
        anl_deriv <- d("a1",p,r1,r0,data,item.focus=item,theta,covariates,cov,num_items,samp_size,num_quadpts)
        p_active <- c(p[paste0("c0_itm",item,"_")],p[grep(paste0("c1_itm",item,"_"),names(p))],p[paste0("a0_itm",item,"_")],p[grep(paste0("a1_itm",item,"_"),names(p))])
        p_current <- p[grep(paste0("a1_itm",item,"_cov",cov),names(p))]
        z <- (anl_deriv[[2]]*p_current - anl_deriv[[1]])/anl_deriv[[2]]
        p_new <- sign(z)*max(abs(z) - tau[t], 0)
        p <- replace(p,names(p_new),p_new)
      }

    #Update and check for convergence: Calculate the difference in parameter estimates from current to previous
    eps = sqrt(sum((p - lastp)^2))

    #Update parameter list
    lastp <- p

    #update the iteration number
    iter = iter + 1
    if(iter == maxit) warning("Coordinate descent iteration limit reached without convergence")

  } #end coordinate descent

 } #end looping through items

  return(p)
}

