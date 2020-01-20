######################
# M-step 1: (Impact) #
######################

Mstep.2pl.impact <-
  function(p,rlist,theta,covariates,maxit,samp_size,num_quadpts){

    p_impact <- c(p[grep("g0",names(p),fixed=T)],p[grep("b0",names(p),fixed=T)])
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
  function(p,rlist,theta,covariates,tau,t,maxit,num_items,samp_size,num_quadpts,num_covariates,anchor){

    #for each item (loop); maximizing i independent logistic regression log-likelihoods (Q functions), with the quadrature points serving as the predictor values
    for (item in 1:num_items) {

      r1 <- rlist[[1]][[item]]
      r0 <- rlist[[2]][[item]]

      anl_deriv <- d("c0",p,r1,r0,data,item.focus=item,theta,covariates,cov=NULL,num_items,samp_size,num_quadpts)
      p_new <- p[paste0("c0_itm",item,"_")] - anl_deriv[[1]]/anl_deriv[[2]]
      p <- replace(p,names(p_new),p_new)

      anl_deriv <- d("a0",p,r1,r0,data,item.focus=item,theta,covariates,cov=NULL,num_items,samp_size,num_quadpts)
      p_new <- p[paste0("a0_itm",item,"_")] - anl_deriv[[1]]/anl_deriv[[2]]
      p <- replace(p,names(p_new),p_new)

      if(any(item == anchor)){
        next
      }

      for(cov in 1:num_covariates){

        if(sum(p[grep(paste0("c1(.*?)cov",cov),names(p))] != 0) > (num_items - 2)){
          next
        }

        anl_deriv <- d("c1",p,r1,r0,data,item.focus=item,theta,covariates,cov,num_items,samp_size,num_quadpts)
        z <- (anl_deriv[[2]]*p[grep(paste0("c1_itm",item,"_cov",cov),names(p),fixed=T)] - anl_deriv[[1]])/anl_deriv[[2]]
        p_new <- sign(z)*max(abs(z) - tau[t], 0)
        p <- replace(p,names(p_new),p_new)
      }

      for(cov in 1:num_covariates){

        if(sum(p[grep(paste0("a1(.*?)cov",cov),names(p))] != 0) > (num_items - 2)){
          next
        }

        anl_deriv <- d("a1",p,r1,r0,data,item.focus=item,theta,covariates,cov,num_items,samp_size,num_quadpts)
        z <- (anl_deriv[[2]]*p[grep(paste0("a1_itm",item,"_cov",cov),names(p),fixed=T)] - anl_deriv[[1]])/anl_deriv[[2]]
        p_new <- sign(z)*max(abs(z) - tau[t], 0)
        p <- replace(p,names(p_new),p_new)
      }

 } #end looping through items

  return(p)
}

