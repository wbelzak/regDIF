#################
# 2PL Traceline #
#################

trace.line.pts <-
  function(p_active,theta,covariates) {
    traceline <- sapply(theta,
                       function(x){
                         1/(1+exp(-((p_active[grep("c0",names(p_active))] + covariates %*% p_active[grep("c1",names(p_active))]) + (p_active[grep("a0",names(p_active))] + covariates %*% p_active[grep("a1",names(p_active))])*x)))})
  }
