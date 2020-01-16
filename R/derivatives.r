#######################
# Partial Derivatives #
#######################


d <-
  function(parm,p,r1,r0,data,item.focus,theta,covariates,cov,num_items,samp_size,num_quadpts) {

    #make space for the trace lines and the derivative E-tables
    logp_deriv <- logq_deriv <- logpq_deriv2 <- matrix(0, nrow = samp_size, ncol = num_quadpts)
    eta.d <- matrix(0, nrow = samp_size, ncol = num_quadpts)

    if(parm == "c0"){
      eta.d <- matrix(1, nrow = samp_size, ncol = num_quadpts)
    } else if(parm == "a0"){
      eta.d <- matrix(rep(matrix(theta, nrow = 1, ncol = length(theta)), samp_size), nrow = samp_size, ncol = length(theta), byrow = TRUE)
    } else if(parm == "c1"){
      eta.d <- matrix(rep(covariates[,cov], num_quadpts), ncol = num_quadpts, nrow = samp_size)
    } else if(parm == "a1"){
      eta.d <- matrix(rep(covariates[,cov], num_quadpts), ncol = num_quadpts, nrow = samp_size)*matrix(rep(matrix(theta, nrow = 1, ncol = length(theta)), samp_size), nrow = samp_size, ncol = length(theta), byrow = TRUE)
    }

    p_active <- c(p[paste0("c0_itm",item.focus,"_")],p[grep(paste0("c1_itm",item.focus,"_"),names(p))],p[paste0("a0_itm",item.focus,"_")],p[grep(paste0("a1_itm",item.focus,"_"),names(p))])
    for(q in 1:num_quadpts){
      logp_deriv[,q] <- (1-trace.line.pts(p_active,theta[q],covariates))*eta.d[,q]
      logq_deriv[,q] <- -trace.line.pts(p_active,theta[q],covariates)*eta.d[,q]
      logpq_deriv2[,q] <- -eta.d[,q]**2*trace.line.pts(p_active,theta[q],covariates)*(1-trace.line.pts(p_active,theta[q],covariates))
    }

    d1 <- sum(logp_deriv*r1) + sum(logq_deriv*r0)
    d2 <- sum(logpq_deriv2*r1) + sum(logpq_deriv2*r0)
    dlist <- list(d1,d2)
  }






