#####################################
# E-step: Evaluating the Q function #
#####################################

Estep.2pl <-
  function(p,data,covariates,theta,t,num_items,samp_size,num_quadpts) { #p is parameters, data is response data, theta is values of latent variable

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


  }


###################
# M-step 2: (DIF) #
###################

Mstep.2pl.dif <-
  function(p,rlist,theta,covariates,tau,t,maxit,num_items,samp_size,num_quadpts,anchor){

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

      fit_dif = optim(par=p_active,fn=ll.2pl.dif,r1=r1,r0=r0,theta=theta,covariates=covariates,tau=tau[t],samp_size=samp_size,num_quadpts=num_quadpts,method="BFGS",control=list(maxit = maxit))

      p <- replace(p,names(fit_dif$par),fit_dif$par)

    }
  return(p)
  }
