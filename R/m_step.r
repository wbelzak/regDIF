##########
# M-step #
##########

Mstep_2pl_dif <-
  function(p,
           responses,
           predictors,
           elist,
           theta,
           itemtypes,
           penalty,
           lambda,
           alpha,
           gamma,
           anchor,
           rasch,
           maxit,
           samp_size,
           num_responses,
           num_items,
           num_quadpts,
           num_predictors) {

  #obtain parameter estimates and posterior probabilities
  p_impact <- c(p[[num_items+1]],p[[num_items+2]])
  etable_all <- elist[[2]]

  #impact mean updates
  for(cov in 1:num_predictors){
    anl_deriv <- d_alpha(p_impact,etable_all,theta,predictors,cov=cov,samp_size,num_items,num_quadpts)
    p_new <- p_impact[grep(paste0("g"),names(p_impact),fixed=T)][cov] - anl_deriv[[1]]/anl_deriv[[2]]
    p_impact <- replace(p_impact,names(p_new),p_new)
  }

  #impact variance updates
  for(cov in 1:num_predictors){
    anl_deriv <- d_phi(p_impact,etable_all,theta,predictors,cov=cov,samp_size,num_items,num_quadpts)
    p_new <- p_impact[grep(paste0("b"),names(p_impact),fixed=T)][cov] - anl_deriv[[1]]/anl_deriv[[2]]
    p_impact <- replace(p_impact,names(p_new),p_new)
  }

  g <- p_impact[grep(paste0("g"),names(p_impact),fixed=T)]
  b <- p_impact[grep(paste0("b"),names(p_impact),fixed=T)]

  p[[num_items+1]] <- replace(p[[num_items+1]],names(g),g)
  p[[num_items+2]] <- replace(p[[num_items+2]],names(b),b)


  #for each item (loop); maximizing i independent logistic regression log-likelihoods (Q functions), with the quadrature points serving as the predictor values
  for (item in 1:num_items) {



    ###################
    # Obtain E-tables #
    ###################

    #get posterior probabilities
    etable <- replicate(n=num_responses[item], elist[[1]][[item]], simplify = F)

    #obtain etables for each response category
    if(itemtypes[item] == "bernoulli" | itemtypes[item] == "categorical"){
      for(resp in 1:num_responses[item]){
        etable[[resp]][which(!(etable[[resp]][,ncol(etable[[resp]])] == resp)),] <- 0
      }
    }

    #get item parameters
    p_item <- p[[item]]
    etable <- lapply(etable, function(x) x[,1:num_quadpts])

    #######################
    # Bernoulli Responses #
    #######################

    if(itemtypes[item] == "bernoulli"){

      #intercept updates
      anl_deriv <- d_bernoulli("c0",p_item,etable,theta,predictors,cov=NULL,samp_size,num_items,num_quadpts)
      p_new <- p_item[grep(paste0("c0_itm",item,"_"),names(p_item),fixed=T)][1] - anl_deriv[[1]]/anl_deriv[[2]]
      p_item <- replace(p_item,names(p_new),p_new)

      #slope updates
      if(rasch == FALSE){
        anl_deriv <- d_bernoulli("a0",p_item,etable,theta,predictors,cov=NULL,samp_size,num_items,num_quadpts)
        p_new <- p_item[grep(paste0("a0_itm",item,"_"),names(p_item),fixed=T)] - anl_deriv[[1]]/anl_deriv[[2]]
        p_item <- replace(p_item,names(p_new),p_new)
      }


      if(!any(item == anchor)){

        p2 <- unlist(p)

        #intercept DIF updates
        for(cov in 1:num_predictors){

          #end routine if only one anchor item is left on each covariate for each item parameter
          if(is.null(anchor) & sum(p2[grep(paste0("c1(.*?)cov",cov),names(p2))] != 0) > (num_items - 2)){
            next
          }

          anl_deriv <- d_bernoulli("c1",p_item,etable,theta,predictors,cov,samp_size,num_items,num_quadpts)
          z <- p_item[grep(paste0("c1_itm",item,"_cov",cov),names(p_item),fixed=T)] - anl_deriv[[1]]/anl_deriv[[2]]
          p_new <- ifelse(penalty == "lasso",soft_threshold(z,alpha,lambda),firm_threshold(z,alpha,lambda,gamma))
          names(p_new) <- names(z)
          p_item <- replace(p_item,names(p_new),p_new)
        }

        #slope DIF updates
        for(cov in 1:num_predictors){

          #end routine if only one anchor item is left on each covariate for each item parameter
          if(is.null(anchor) & sum(p2[grep(paste0("a1(.*?)cov",cov),names(p2))] != 0) > (num_items - 2)){
            next
          }

          if(rasch == FALSE){
            anl_deriv <- d_bernoulli("a1",p_item,etable,theta,predictors,cov,samp_size,num_items,num_quadpts)
            z <- p_item[grep(paste0("a1_itm",item,"_cov",cov),names(p_item),fixed=T)] - anl_deriv[[1]]/anl_deriv[[2]]
            p_new <- ifelse(penalty == "lasso",soft_threshold(z,alpha,lambda),firm_threshold(z,alpha,lambda,gamma))
            names(p_new) <- names(z)
            p_item <- replace(p_item,names(p_new),p_new)
          }
        }
      }

      #########################
      # Categorical Responses #
      #########################

  } else if(itemtypes[item] == "categorical"){

      #intercept updates
      anl_deriv <- d_categorical("c0",p_item,etable,theta,predictors,thr=NULL,cov=NULL,samp_size,num_responses[[item]],num_items,num_quadpts)
      p_new <- p_item[grep(paste0("c0_itm",item,"_"),names(p_item),fixed=T)][1] - anl_deriv[[1]]/anl_deriv[[2]]
      p_item <- replace(p_item,names(p_new),p_new)

      #threshold updates
      if(num_responses[item] > 2){
        for(thr in 2:(num_responses[item]-1)){
          anl_deriv <- d_categorical("c0",p_item,etable,theta,predictors,thr=thr,cov=NULL,samp_size,num_responses[[item]],num_items,num_quadpts)
          p_new <- p_item[grep(paste0("c0_itm",item,"_"),names(p_item),fixed=T)][thr] - anl_deriv[[1]]/anl_deriv[[2]]
          p_item <- replace(p_item,names(p_new),p_new)
        }
      }

      #slope updates
      if(rasch == FALSE){ #skip slope estimate updates if rasch is TRUE
        anl_deriv <- d_categorical("a0",p_item,etable,theta,predictors,thr=NULL,cov=NULL,samp_size,num_responses[[item]],num_items,num_quadpts)
        p_new <- p_item[grep(paste0("a0_itm",item,"_"),names(p_item),fixed=T)] - anl_deriv[[1]]/anl_deriv[[2]]
        p_item <- replace(p_item,names(p_new),p_new)
      }

      if(!any(item == anchor)){

        p2 <- unlist(p) #unlist to check for anchor identification

        #intercept DIF updates
        for(cov in 1:num_predictors){

          #end routine if only one anchor item is left on each covariate for each item parameter
          if(is.null(anchor) & sum(p2[grep(paste0("c1(.*?)cov",cov),names(p2))] != 0) > (num_items - 2)){
            next
          }

          anl_deriv <- d_categorical("c1",p_item,etable,theta,predictors,thr=NULL,cov,samp_size,num_responses[[item]],num_items,num_quadpts)
          z <- p_item[grep(paste0("c1_itm",item,"_cov",cov),names(p_item),fixed=T)] - anl_deriv[[1]]/anl_deriv[[2]]
          p_new <- ifelse(penalty == "lasso",soft_threshold(z,alpha,lambda),firm_threshold(z,alpha,lambda,gamma))
          names(p_new) <- names(z)
          p_item <- replace(p_item,names(p_new),p_new)
        }

        #slope DIF updates
        for(cov in 1:num_predictors){

          #end routine if only one anchor item is left on each covariate for each item parameter
          if(is.null(anchor) & sum(p2[grep(paste0("a1(.*?)cov",cov),names(p2))] != 0) > (num_items - 2)){
            next
          }

          if(rasch == FALSE){
            anl_deriv <- d_categorical("a1",p_item,etable,theta,predictors,thr=NULL,cov,samp_size,num_responses[[item]],num_items,num_quadpts)
            z <- p_item[grep(paste0("a1_itm",item,"_cov",cov),names(p_item),fixed=T)] - anl_deriv[[1]]/anl_deriv[[2]]
            p_new <- ifelse(penalty == "lasso",soft_threshold(z,alpha,lambda),firm_threshold(z,alpha,lambda,gamma))
            names(p_new) <- names(z)
            p_item <- replace(p_item,names(p_new),p_new)
          }
        }
      }


      ######################
      # Gaussian Responses #
      ######################

    } else if(itemtypes[item] == "gaussian"){

      #intercept updates
      anl_deriv <- d_mu_gaussian("c0",p_item,etable,theta,responses[,item],predictors,cov=NULL,samp_size,num_items,num_quadpts)
      p_new <- p_item[grep(paste0("c0_itm",item,"_"),names(p_item),fixed=T)][1] - anl_deriv[[1]]/anl_deriv[[2]]
      p_item <- replace(p_item,names(p_new),p_new)

      #slope updates
      if(rasch == FALSE){ #skip slope estimate updates if rasch is TRUE
        anl_deriv <- d_mu_gaussian("a0",p_item,etable,theta,responses[,item],predictors,cov=NULL,samp_size,num_items,num_quadpts)
        p_new <- p_item[grep(paste0("a0_itm",item,"_"),names(p_item),fixed=T)] - anl_deriv[[1]]/anl_deriv[[2]]
        p_item <- replace(p_item,names(p_new),p_new)
      }

      #residual updates
      anl_deriv <- d_sigma_gaussian("s0",p_item,etable,theta,responses[,item],predictors,cov=NULL,samp_size,num_items,num_quadpts)
      p_new <- p_item[grep(paste0("s_itm",item,"_"),names(p_item),fixed=T)][1] - anl_deriv[[1]]/anl_deriv[[2]]
      p_item <- replace(p_item,names(p_new),p_new)


      if(!any(item == anchor)){

        #residual dif updates
        for(cov in 1:num_predictors){
          anl_deriv <- d_sigma_gaussian("s1",p_item,etable,theta,responses[,item],predictors,cov=cov,samp_size,num_items,num_quadpts)
          p_new <- p_item[grep(paste0("s_itm",item,"_cov",cov),names(p_item),fixed=T)][1] - anl_deriv[[1]]/anl_deriv[[2]]
          p_item <- replace(p_item,names(p_new),p_new)
        }

        p2 <- unlist(p)

        #intercept DIF updates
        for(cov in 1:num_predictors){

          #end routine if only one anchor item is left on each covariate for each item parameter
          if(is.null(anchor) & sum(p2[grep(paste0("c1(.*?)cov",cov),names(p2))] != 0) > (num_items - 2)){
            next
          }

          anl_deriv <- d_mu_gaussian("c1",p_item,etable,theta,responses[,item],predictors,cov,samp_size,num_items,num_quadpts)
          z <- p_item[grep(paste0("c1_itm",item,"_cov",cov),names(p_item),fixed=T)] - anl_deriv[[1]]/anl_deriv[[2]]
          p_new <- ifelse(penalty == "lasso",soft_threshold(z,alpha,lambda),firm_threshold(z,alpha,lambda,gamma))
          names(p_new) <- names(z)
          p_item <- replace(p_item,names(p_new),p_new)
        }

        #slope DIF updates
        for(cov in 1:num_predictors){

          #end routine if only one anchor item is left on each covariate for each item parameter
          if(is.null(anchor) & sum(p2[grep(paste0("a1(.*?)cov",cov),names(p2))] != 0) > (num_items - 2)){
            next
          }

          if(rasch == FALSE){
            anl_deriv <- d_mu_gaussian("a1",p_item,etable,theta,responses[,item],predictors,cov,samp_size,num_items,num_quadpts)
            z <- p_item[grep(paste0("a1_itm",item,"_cov",cov),names(p_item),fixed=T)] - anl_deriv[[1]]/anl_deriv[[2]]
            p_new <- ifelse(penalty == "lasso",soft_threshold(z,alpha,lambda),firm_threshold(z,alpha,lambda,gamma))
            names(p_new) <- names(z)
            p_item <- replace(p_item,names(p_new),p_new)
          }
        }
      }

    }


    p[[item]] <- replace(p[[item]],names(p_item),p_item)


  } #end looping through items

  return(p)

}

