#################
# 2PL Traceline #
#################

trace.line.pts <-
  function(p_active,theta,covariates) {
    active_cov_c <- as.matrix(covariates[,as.numeric(gsub(".*([0-9]+)", "\\1", names(p_active[grep("c1",names(p_active))])))])
    active_cov_a <- as.matrix(covariates[,as.numeric(gsub(".*([0-9]+)", "\\1", names(p_active[grep("a1",names(p_active))])))])

    traceline <- 1/(1+exp(-((p_active[grep("c0",names(p_active))] + active_cov_c %*% p_active[grep("c1",names(p_active))]) + (p_active[grep("a0",names(p_active))] + active_cov_a %*% p_active[grep("a1",names(p_active))])*theta))) #logistic CDF

  }
