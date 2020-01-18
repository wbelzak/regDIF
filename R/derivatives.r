#######################
# Partial Derivatives #
#######################

d <-
  function(parm,p,r1,r0,data,item.focus,theta,covariates,cov,num_items,samp_size,num_quadpts) {

    #make space for the trace lines and the derivative E-tables
    traceline <- matrix(0, nrow = samp_size, ncol = num_quadpts)
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
    traceline <- trace.line.pts(p_active,theta,covariates)

    d1 <- sum((1-traceline)*eta.d*r1) + sum(-traceline*eta.d*r0)
    d2 <- sum(-traceline*(1-traceline)*eta.d**2*r1) + sum(-traceline*(1-traceline)*eta.d**2*r0)
    dlist <- list(d1,d2)
  }




