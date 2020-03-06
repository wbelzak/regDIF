######################
# M-step 1: (Impact) #
######################

Mstep_2pl_impact <-
  function(p,
           elist,
           theta,
           predictors,
           maxit,
           samp_size,
           num_items){

    #obtain parameter estimates and posterior probabilities
    p_impact <- c(p[[num_items+1]],p[[num_items+2]])
    etable_all <- elist[[2]]

    #update impact parameters
    fit_impact = optim(par=p_impact,fn=ll.2pl.impact,etable_all=etable_all,theta=theta,predictors=predictors,samp_size=samp_size,method="BFGS",control=list(maxit = maxit))
    g <- fit_impact$par[grep("g",names(fit_impact$par),fixed=T)]
    b <- fit_impact$par[grep("b",names(fit_impact$par),fixed=T)]

    #collect impact parameter estimates
    p[[num_items+1]] <- replace(p[[num_items+1]],names(g),g)
    p[[num_items+2]] <- replace(p[[num_items+2]],names(b),b)

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
           itemtypes,
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

      ###################
      # Obtain E-tables #
      ###################

      #get posterior probabilities
      etable <- replicate(n=num_responses[item], elist[[1]][[item]], simplify = F)

      #obtain etables for each response category
      if(itemtypes[item] == "categorical"){
        for(resp in 1:num_responses[item]){
          etable[[resp]][which(!(etable[[resp]][,ncol(etable[[resp]])] == resp)),] <- 0
        }
      }

      #get item parameters
      p_item <- p[[item]]
      etable <- lapply(etable, function(x) x[,1:num_quadpts])

      #################################
      # Mu (mean) baseline parameters #
      #################################

      #intercept updates
      anl_deriv <- d_mu("c0",p_item,etable,theta,responses[,item],predictors,itemtypes[item],thr=NULL,cov=NULL,samp_size,num_responses[[item]],num_items,num_quadpts)
      p_new <- p_item[grep(paste0("c0_itm",item,"_"),names(p_item),fixed=T)][1] - anl_deriv[[1]]/anl_deriv[[2]]
      p_item <- replace(p_item,names(p_new),p_new)

      #threshold updates for graded response model
      if(num_responses[item] > 2){
        for(thr in 2:(num_responses[item]-1)){
          anl_deriv <- d_mu("c0",p_item,etable,theta,responses[,item],predictors,itemtypes[item],thr,cov=NULL,samp_size,num_responses[[item]],num_items,num_quadpts)
          p_new <- p_item[grep(paste0("c0_itm",item,"_"),names(p_item),fixed=T)][thr] - anl_deriv[[1]]/anl_deriv[[2]]
          p_item <- replace(p_item,names(p_new),p_new)
        }
      }

      #slope updates
      if(rasch == FALSE){ #skip slope estimate updates if rasch is TRUE
        anl_deriv <- d_mu("a0",p_item,etable,theta,responses[,item],predictors,itemtypes[item],thr=NULL,cov=NULL,samp_size,num_responses[[item]],num_items,num_quadpts)
        p_new <- p_item[grep(paste0("a0_itm",item,"_"),names(p_item),fixed=T)] - anl_deriv[[1]]/anl_deriv[[2]]
        p_item <- replace(p_item,names(p_new),p_new)
      }

      #########################################
      # Residual variance baseline parameters #
      #########################################

      if(itemtypes[item] == "continuous"){

        anl_deriv <- d_sigma("s0",p_item,etable,theta,responses[,item],predictors,itemtypes[item],cov=NULL,samp_size,num_items,num_quadpts)
        p_new <- p_item[grep(paste0("s_itm",item,"_"),names(p_item),fixed=T)][1] - anl_deriv[[1]]/anl_deriv[[2]]
        # p_new <- ifelse(p_new <= 0, 0.01, p_new)
        p_item <- replace(p_item,names(p_new),p_new)

      }


      #skip DIF estimate updates if anchor item
      if(any(item == anchor)){
        p[[item]] <- replace(p[[item]],names(p_item),p_item)
        next
      }

      ############################
      # Mu (mean) DIF parameters #
      ############################

      p2 <- unlist(p) #unlist to check for anchor identification

      #intercept DIF updates
      for(cov in 1:num_predictors){

        #end routine if only one anchor item is left on each covariate for each item parameter
        if(is.null(anchor) & sum(p2[grep(paste0("c1(.*?)cov",cov),names(p2))] != 0) > (num_items - 2)){
          next
        }

        anl_deriv <- d_mu("c1",p_item,etable,theta,responses[,item],predictors,itemtypes[item],thr=NULL,cov,samp_size,num_responses[[item]],num_items,num_quadpts)
        z <- (anl_deriv[[2]]*p_item[grep(paste0("c1_itm",item,"_cov",cov),names(p_item),fixed=T)] - anl_deriv[[1]])/anl_deriv[[2]]
        p_new <- sign(z)*max(abs(z) - penalty[pen], 0)
        p_item <- replace(p_item,names(p_new),p_new)
      }

      #slope DIF updates
      for(cov in 1:num_predictors){

        #end routine if only one anchor item is left on each covariate for each item parameter
        if(is.null(anchor) & sum(p2[grep(paste0("a1(.*?)cov",cov),names(p2))] != 0) > (num_items - 2)){
          next
        }

        if(rasch == FALSE){
          anl_deriv <- d_mu("a1",p_item,etable,theta,responses[,item],predictors,itemtypes[item],thr=NULL,cov,samp_size,num_responses[[item]],num_items,num_quadpts)
          z <- (anl_deriv[[2]]*p_item[grep(paste0("a1_itm",item,"_cov",cov),names(p_item),fixed=T)] - anl_deriv[[1]])/anl_deriv[[2]]
          p_new <- sign(z)*max(abs(z) - penalty[pen], 0)
          p_item <- replace(p_item,names(p_new),p_new)
        }
      }

      ####################################
      # Residual variance DIF parameters #
      ####################################

      if(itemtypes[item] == "continuous"){

        for(cov in 1:num_predictors){
          # anl_deriv <- d_sigma("s1",p_item,etable,theta,responses[,item],predictors,itemtypes,cov=cov,samp_size,num_items,num_quadpts)
          # z <- (anl_deriv[[2]]*p_item[grep(paste0("s_itm",item,"_cov",cov),names(p_item),fixed=T)] - anl_deriv[[1]])/anl_deriv[[2]]
          # p_new <- sign(z)*max(abs(z) - penalty[pen], 0)
          # p_item <- replace(p_item,names(p_new),p_new)

          anl_deriv <- d_sigma("s1",p_item,etable,theta,responses[,item],predictors,itemtypes[item],cov=cov,samp_size,num_items,num_quadpts)
          p_new <- p_item[grep(paste0("s_itm",item,"_cov",cov),names(p_item),fixed=T)][1] - anl_deriv[[1]]/anl_deriv[[2]]
          p_item <- replace(p_item,names(p_new),p_new)
        }

      }


      #collect all updated parameters
      p[[item]] <- replace(p[[item]],names(p_item),p_item)

 } #end looping through items

  return(p)
}

