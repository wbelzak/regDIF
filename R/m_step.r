######################
# M-step 1: (Impact) #
######################

Mstep_2pl_impact <-
  function(p,
           elist,
           theta,
           predictors,
           maxit,
           samp_size){

    #obtain parameter estimates and posterior probabilities
    p_impact <- c(p[[7]],p[[8]])
    etable_all <- elist[[2]]

    #M-step 1: (Impact)
    fit_impact = optim(par=p_impact,fn=ll.2pl.impact,etable_all=etable_all,theta=theta,predictors=predictors,samp_size=samp_size,method="BFGS",control=list(maxit = maxit))
    g <- fit_impact$par[grep("g",names(fit_impact$par),fixed=T)]
    b <- fit_impact$par[grep("b",names(fit_impact$par),fixed=T)]

    p[[7]] <- replace(p[[7]],names(g),g)
    p[[8]] <- replace(p[[8]],names(b),b)

    ll <- fit_impact$value

    return(list(p,ll))

  }

###################
# M-step 2: (DIF) #
###################
Mstep_2pl_dif <-
  function(p,
           responses,
           predictors,
           elist,
           theta,
           penalty,
           pen,
           anchor,
           rasch,
           maxit,
           samp_size,
           num_responses,
           num_items,
           num_quadpts,
           num_predictors){

    #for each item (loop); maximizing i independent logistic regression log-likelihoods (Q functions), with the quadrature points serving as the predictor values
    for (item in 1:num_items) {

      #get etable information
      etable <- replicate(n=num_responses[item], elist[[1]][[item]], simplify = F)
      for(resp in 1:num_responses[item]){
        etable[[resp]][which(!(etable[[resp]][,ncol(etable[[resp]])] == resp)),] <- 0
      }
      p_item <- p[[item]]
      etable <- lapply(etable, function(x) x[,1:num_quadpts])

      anl_deriv <- d("c0",p_item,etable,theta,predictors,thr=NULL,cov=NULL,samp_size,num_responses[[item]],num_items,num_quadpts)
      p_new <- p_item[grep(paste0("c0_itm",item,"_"),names(p_item),fixed=T)][1] - anl_deriv[[1]]/anl_deriv[[2]]
      p_item <- replace(p_item,names(p_new),p_new)

      if(num_responses[item] > 2){
        for(thr in 2:(num_responses[item]-1)){
          anl_deriv <- d("c0",p_item,etable,theta,predictors,thr,cov=NULL,samp_size,num_responses[[item]],num_items,num_quadpts)
          p_new <- p_item[grep(paste0("c0_itm",item,"_"),names(p_item),fixed=T)][thr] - anl_deriv[[1]]/anl_deriv[[2]]
          p_item <- replace(p_item,names(p_new),p_new)
        }
      }

      if(rasch == FALSE){
        anl_deriv <- d("a0",p_item,etable,theta,predictors,thr=NULL,cov=NULL,samp_size,num_responses[[item]],num_items,num_quadpts)
        p_new <- p_item[grep(paste0("a0_itm",item,"_"),names(p_item),fixed=T)] - anl_deriv[[1]]/anl_deriv[[2]]
        p_item <- replace(p_item,names(p_new),p_new)
      }

      if(any(item == anchor)){
        p[[item]] <- replace(p[[item]],names(p_item),p_item)
        next
      }

      p2 <- unlist(p)
      for(cov in 1:num_predictors){

        #end routine if only one anchor item is left on each covariate for each item parameter
        if(is.null(anchor) & sum(p2[grep(paste0("c1(.*?)cov",cov),names(p2))] != 0) > (num_items - 2)){
          next
        }

        anl_deriv <- d("c1",p_item,etable,theta,predictors,thr=NULL,cov,samp_size,num_responses[[item]],num_items,num_quadpts)
        z <- (anl_deriv[[2]]*p_item[grep(paste0("c1_itm",item,"_cov",cov),names(p_item),fixed=T)] - anl_deriv[[1]])/anl_deriv[[2]]
        p_new <- sign(z)*max(abs(z) - penalty[pen], 0)
        p_item <- replace(p_item,names(p_new),p_new)
      }

      for(cov in 1:num_predictors){

        #end routine if only one anchor item is left on each covariate for each item parameter
        if(is.null(anchor) & sum(p2[grep(paste0("a1(.*?)cov",cov),names(p2))] != 0) > (num_items - 2)){
          next
        }

        anl_deriv <- d("a1",p_item,etable,theta,predictors,thr=NULL,cov,samp_size,num_responses[[item]],num_items,num_quadpts)
        z <- (anl_deriv[[2]]*p_item[grep(paste0("a1_itm",item,"_cov",cov),names(p_item),fixed=T)] - anl_deriv[[1]])/anl_deriv[[2]]
        p_new <- sign(z)*max(abs(z) - penalty[pen], 0)
        p_item <- replace(p_item,names(p_new),p_new)
      }

      p[[item]] <- replace(p[[item]],names(p_item),p_item)
 } #end looping through items

  return(p)
}

