#' Regularized Differential Item Functioning
#'
#' Regularization of DIF parameters in item response theory (IRT) and moderated nonlinear factor analysis (MNLFA) models.
#'
#' @param y Matrix or dataframe of item responses. Supports dichotomous (e.g., 0,1), ordinal (e.g., 0,1,2,...,m), and continuous item responses. See more information in \code{itemtypes} description.
#' @param x Matrix or dataframe of explanatory predictors (i.e., DIF covariates). Supports categorical and continuous predictors.
#' @param itemtypes Character value or vector indicating the item response distributions. For scales where item responses are of one type only, the user may input one character value indicating the type (e.g., \code{itemtypes = "categorical"}). For mixed item types, the user must specify a vector of characters in the order that corresponds to the response matrix (e.g. \code{itemtypes = c(rep("categorical",3),rep("continuous",3))}). Supports categorical and continuous responses, including Bernoulli with logistic link (i.e., binary outcomes), categorical with ordered logistic link (i.e., graded response model), and Gaussian with identity link (i.e., factor analysis).
#' @param penalty Numeric vector of tuning parameter values \eqn{\ge} 0. Must be in descending order, from largest to smallest values (e.g., \code{seq(1,0,-.01)}.
#' @param anchor Optional numeric vector indicating which items are anchors (e.g., \code{anchor = 1}). Default is \code{NULL}, meaning at least one DIF parameter (per covariate) will be fixed to zero as the penalty values approach 0. This is required to identify the model.
#' @param standardize Logical value indicating whether to normalize the predictors. Default is \code{TRUE}.
#' @param rasch Logical value indicating whether to constrain item slopes to 1 (i.e., equal slopes). If \code{TRUE}, no slope DIF will be evaluated. Default is \code{FALSE}.
#'
#' @return Function returns an object of class \code{regDIF}
#'
#' @examples
#' \dontrun{
#'
#' library(regDIF)
#' head(ida)
#' y <- ida[,1:6]
#' x <- ida[,7:9]
#' penalty <- seq(1,0,-.05)
#' fit <- regDIF(y, x, penalty)
#' fit
#'
#' }
#'
#' @import stats utils
#'
#' @export

regDIF <- function(y,
                   x,
                   itemtypes = c("categorical","continuous"),
                   penalty,
                   anchor = NULL,
                   standardize = TRUE,
                   rasch = FALSE){

  ##############
  # Preprocess #
  ##############
  # library(lavaan)
  # y <- HolzingerSwineford1939[,7:15]
  # x <- HolzingerSwineford1939[,2]
  # y[,1] <- ifelse(y[,1]>5,2,1)
  # y[,2] <- as.numeric(cut(y[,2],seq(min(y[,2]),max(y[,2]),by=1.75)))
  # write.table(cbind(1:nrow(y),y,x),"C:\\Users\\wbelz\\Dropbox\\Will\\Research\\Regularized IRT\\mixed_data2.dat",col.names = F,row.names = F)
  # itemtypes = c(rep('categorical',2),rep('continuous',7))
  # penalty = 0
  # anchor = 1
  # standardize = F
  # num_quadpts=51
  # rasch =F


  #item response warnings
  # if(responsetype == "binary" | any(responsetype == "binary")) {if(any(!(responses == 0 | responses == 1), na.rm = TRUE)) stop("Item responses must be scored 1 for yes/correct and 0 for no/incorrect.", call. = TRUE)}
  if(!(any(itemtypes == "categorical") | any(itemtypes == "continuous"))) stop("Item types must be either 'categorical' or 'continuous'.")
  #penalty warnings
  if(any(penalty < 0)) stop("Penalty values must be non-negative.", call. = TRUE)
  if(length(penalty) > 1 & all(diff(penalty) >= 0)) stop("Penalty values must be in descending order (e.g., penalty = c(-2,-1,0)).")
  #anchor warnings
  if(is.null(anchor) & length(penalty) == 1){if(penalty == 0) stop("Anchor item must be specified without penalty.", call. = TRUE)}
  if(!is.null(anchor) & !is.numeric(anchor)) stop("Anchor items must be numeric (e.g., anchor = 1).", call. = TRUE)
  #quadrature warnings
  # if(quadpts < 1) stop("The number of quadrature points must be greater than 1.\n  15 or more points may be required to get accurate estimates.")
  # if(!is.numeric(quadpts)) stop("Quadrature points must be numeric (e.g., quadpts = 15).")

  #data
  responses <- y
  predictors <- x

  #get number of quadrature points
  num_quadpts <- 51

  #get latent variable values (i.e., predictor values) for quadrature and tracelines
  theta <- seq(-6, 6, length.out = num_quadpts)


  #turn data into numeric if not already
  if(any(!sapply(responses,function(x)is.numeric(x)))){responses <- sapply(responses,function(x)as.numeric(x))}
  if(any(!sapply(predictors,function(x)is.numeric(x)))){predictors <- sapply(predictors,function(x)as.numeric(x))}


  if(length(itemtypes) == 1){
    itemtypes <- rep(itemtypes, num_items)
  }

  if(any(itemtypes == "categorical")){
    responses[,which(itemtypes == "categorical")] <- sapply(responses[,which(itemtypes == "categorical")], function(x) as.numeric(as.factor(x)))
  }

  #speeds up computation
  responses <- as.matrix(responses)
  predictors <- as.matrix(predictors)
  num_items <- dim(responses)[2]
  samp_size <- dim(responses)[1]
  num_predictors <- dim(predictors)[2]

  #get item response types
  num_responses <- rep(1,num_items)
  if(any(itemtypes == "categorical")){
    num_responses[which(itemtypes == "categorical")] <- apply(responses[,which(itemtypes == "categorical")], 2, function(x) length(unique(na.omit(x))))
  }

  #standardize predictors
  if(standardize == TRUE){
    predictors <- scale(predictors)
  }

  ###################
  # Starting Values #
  ###################
  p <- replicate(n=num_items+2,list(NA),simplify=F)
  for(item in 1:num_items){
    if(num_responses[item] > 2){ #categorical item responses
      p[[item]] <- c(0, seq(.25,1,length.out = num_responses[item]-2), 1, rep(0, num_predictors), rep(0, num_predictors))
      names(p[[item]]) <- c(paste0('c0_itm',item,"_int"),paste0('c0_itm',item,"_thr",1:(num_responses[item]-2),"_"),paste0('a0_itm',item,"_"),paste0('c1_itm',item,"_cov",1:num_predictors),paste0('a1_itm',item,"_cov",1:num_predictors))
    } else if(num_responses[item] == 2){ #bernoulli item responses
      p[[item]] <- c(0, 1, rep(0, num_predictors), rep(0, num_predictors))
      names(p[[item]]) <- c(paste0('c0_itm',item,"_int"),paste0('a0_itm',item,"_"),paste0('c1_itm',item,"_cov",1:num_predictors),paste0('a1_itm',item,"_cov",1:num_predictors))
    } else if(num_responses[item] == 1){ #normal item responses
      p[[item]] <- c(mean(responses[,item]), sqrt(.5*var(responses[,item])), rep(0, num_predictors), rep(0, num_predictors), .5*var(responses[,item]), rep(0, num_predictors))
      names(p[[item]]) <- c(paste0('c0_itm',item,"_int"),paste0('a0_itm',item,"_"),paste0('c1_itm',item,"_cov",1:num_predictors),paste0('a1_itm',item,"_cov",1:num_predictors),paste0('s_itm',item,"_"),paste0('s_itm',item,"_cov",1:num_predictors))
    }
  }
  p[[(num_items+1)]] <- p[[(num_items+2)]] <- rep(0,num_predictors)
  names(p[[(num_items+1)]]) <- paste0(rep(paste0('g',1:num_predictors)))
  names(p[[(num_items+2)]]) <- paste0(rep(paste0('b',1:num_predictors)))
  call <- match.call()
  final <- list(Penalty = rep(NA,length(penalty)),
                BIC = rep(NA,length(penalty)),
                Impact = replicate(n=length(penalty), NA, simplify = F),
                DIF = replicate(n=length(penalty), NA, simplify = F),
                call = call)



  ##################################
  # Reg-DIF - Loop through penalty #
  ##################################
  message(paste0('Running ',length(penalty),if(length(penalty) == 1){' model.'}else{' models.'}))

  for(pen in 1:length(penalty)){

    #Maximization settings
    lastp <- p
    eps <- Inf
    tol = 10^-5
    maxit = 2000
    iter = 0

    ################
    # EM Algorithm #
    ################

    #loop until convergence or maximum number of iterations
    while(eps > tol & iter < maxit){

      #E-step: Evaluate Q function with current parameter estimates p
      elist <- Estep_2pl(p,responses,predictors,theta,itemtypes,samp_size,num_items,num_responses,num_quadpts)

      #M-step 1: Optimize impact parameters
      lv <- Mstep_2pl_impact(p,elist,theta,predictors,maxit,samp_size,num_items)

      #M-step 2: Optimize DIF parameters
      p <- Mstep_2pl_dif(lv[[1]],responses,predictors,elist,theta,penalty,itemtypes,pen,anchor,rasch,maxit,samp_size,num_responses,num_items,num_quadpts,num_predictors)

      #Update and check for convergence: Calculate the difference in parameter estimates from current to previous
      eps = sqrt(sum((unlist(p)-unlist(lastp))^2))

      #Update parameter list
      lastp <- p

      #update the iteration number
      iter = iter + 1
      if(iter == maxit) warning("EM iteration limit reached without convergence")

      cat('\r',"Models Completed:",pen-1," Iteration:",iter," Change:",round(eps,7)) #print information about optimization
      flush.console()

      # print(p[[2]])
    } #End EM once converged or reached iteration limit

    ###############
    # Postprocess #
    ###############

    #obtain likelihood values for all items (** Currently works for only binary items **)
    itemtrace <- lapply(num_responses, function(x) replicate(n=x, matrix(0,nrow=samp_size,ncol=num_quadpts), simplify = F))
    ll_dif <- ll_dif_pen <- current_pen <- replicate(n=num_items, 0, simplify = F)
    for (item in 1:num_items) { #loop through items
      etable1 <- etable2 <- elist[[1]][[item]]
      etable1[which(!(etable1[,ncol(etable1)] == 1)),] <- 0
      etable2[which(!(etable2[,ncol(etable2)] == 2)),] <- 0
      p_item <- p[[item]]
      etable1 <- etable1[,1:num_quadpts]
      etable2 <- etable2[,1:num_quadpts]

      itemtrace[[item]] <- categorical_traceline_pts(p[[item]],theta,predictors,samp_size,num_responses[item],num_quadpts) #computing probability of endorsement for theta value using current estimate of a and b
      ll_dif[[item]] <- (-1)*(sum(etable2*log(itemtrace[[item]][[2]]), na.rm = TRUE) + sum(etable1*log(itemtrace[[item]][[1]]), na.rm = TRUE))
      current_pen <- (penalty[pen]*sum(c(abs(p_item[grep("c1_itm",names(p_item))]), abs(p_item[grep("a1_itm",names(p_item))])), na.rm = TRUE))
      ll_dif_pen[[item]] <- ll_dif[[item]] - current_pen #Q function we want to minimize
    }

    p2 <- unlist(p)
    ll <- lv[[2]] + sum(unlist(ll_dif_pen))
    bic <- length(p2[p2!=0])*log(samp_size) + 2*ll


    #Organize parameters into presentable form
    parms_impact <- rbind(p[[num_items+1]],p[[num_items+2]])
    if(is.null(colnames(x)) | length(colnames(x)) == 0){
      impact_colnames <- c(paste0('cov',1:num_predictors))
    } else{
      impact_colnames <- colnames(x)
    }
    colnames(parms_impact) <- impact_colnames
    rownames(parms_impact) <- c('Mean','Variance')


    parms_dif <- cbind(p2[grep("c0_itm",names(p2))],
                       do.call(rbind, split(p2[grep("c1_itm",names(p2),fixed=T)],rep(1:num_items,each=num_predictors))),
                       p2[grep("a0_itm",names(p2))],
                       do.call(rbind, split(p2[grep("a1_itm",names(p2),fixed=T)],rep(1:num_items,each=num_predictors))))
    #get DIF column names
    if(is.null(colnames(x)) | length(colnames(x)) == 0){
      dif_colnames_int <- paste0('int_cov',1:num_predictors)
      dif_colnames_slp <- paste0('slp_cov',1:num_predictors)
    } else{
      dif_colnames_int <- paste0('int_',colnames(x))
      dif_colnames_slp <- paste0('slp_',colnames(x))
    }
    #get DIF row names
    if(is.null(rownames(y)) | length(rownames(y)) == 0){
      dif_rownames <- paste0(1:num_items)
    } else{
      dif_rownames <- colnames(y)
    }

    colnames(parms_dif) <- c('int_base',
                             dif_colnames_int,
                             'slp_base',
                             dif_colnames_slp)
    rownames(parms_dif) <- dif_rownames

    final$Penalty[pen] <- penalty[pen]
    final$BIC[pen] <- round(bic,2)
    final$Impact[[pen]] <- round(parms_impact,2)
    final$DIF[[pen]] <- round(parms_dif,2)


    #stop if penalty is too small on first run (this leads to different results b/c of identification constraints on DIF parameters)
    if(is.null(anchor) & pen == 1 & sum(abs(p2[grep(paste0("cov"),names(p2))])) > 0){
      print(coef(final))
      stop("First penalty value is too small.\n  Two Options:\n  1. Increase first penalty value large enough to ensure all DIF parameters are removed from the model.\n  2. Provide anchor item(s).", call. = TRUE)
    }

    #stop if there is a large change in DIF parameters
    if(pen > 1){
      second_last <- sum(final$DIF[[pen-1]][,-grep("0",colnames(final$DIF[[pen-1]]))] == 0)
      last <- sum(final$DIF[[pen]][,-grep("0",colnames(final$DIF[[pen-1]]))] == 0)
      if((second_last - last) > (num_predictors*num_items)){
        print(final)
        stop(paste0("Large increase in the number of DIF parameters from iteration ",pen-1," to ",pen,".\n  Two Options:\n  1. Provide smaller differences between penalty values.\n  2. Provide anchor item(s)."), call. = TRUE)
      }
    }

    cat('\r',"Models Completed:",pen," Iteration:",iter," Change:",round(eps,7)) #print information about optimization
    flush.console()

  } #Terminate Reg-DIF



#Obtain final results
message(paste0('\nFinished.'))
class(final) <- "regDIF"
return(final)

}

