#' Regularized Differential Item Functioning
#'
#' Regularization of DIF parameters in item response theory (IRT) and moderated nonlinear factor analysis (MNLFA) models.
#'
#' @param x Matrix or dataframe of predictors (i.e., DIF covariates). Supports categorical and continuous predictors.
#' @param y Matrix or dataframe of item responses. Supports Bernoulli (e.g., 0,1), categorical (e.g., 0,1,2,...,k), and Gaussian item responses.
#' @param itemtypes Character value or vector indicating the item response distributions. For scales where item responses are of one type only, the user may input one character value indicating the type (e.g., \code{"categorical"}). For mixed item types, the user must specify a vector of characters in the order that corresponds to the response matrix (e.g. \code{c(rep("categorical",3),rep("gaussian",3))}). Supports categorical and continuous responses, including Bernoulli with logistic link (i.e., binary outcomes), categorical with ordered logistic link (i.e., graded response model), and Gaussian with identity link (i.e., factor analysis).
#' @param penalty Character value indicating the penalty function to use. Supports the least absolute selection and shrinkage operator (lasso) and the minimax concave penalty (mcp).
#' @param nlambda Numeric value indicating the number of lambda values to fit.
#' @param lambda.max Numberic value indicating the maximum lambda (shrinkage tuning parameter) to fit. Default is 3.
#' @param gamma Numeric value indicating the gamma (tapering tuning parameter) to fit. Default is 3.
#' @param lambda Optional numeric vector of tuning parameter values \eqn{\ge} 0. If supplied, must be in descending order, from largest to smallest values (e.g., \code{seq(1,0,-.01)}.
#' @param anchor Optional numeric vector indicating which items are anchors (e.g., \code{anchor = 1}). Default is \code{NULL}, meaning at least one DIF parameter (per covariate) will be fixed to zero as lambda approaches 0. This is required to identify the model.
#' @param rasch Logical value indicating whether to constrain item slopes to 1 (i.e., equal slopes). If \code{TRUE}, no slope DIF will be evaluated. Default is \code{FALSE}.
#' @param control Optional list of control parameters. See documentation.
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
#' fit <- regDIF(x, y, itemtypes = "bernoulli", penalty = "lasso")
#' fit
#'
#' }
#'
#' @import stats utils
#' @importFrom Rcpp sourceCpp
#' @useDynLib regDIF, .registration = TRUE
#'
#' @export

regDIF <- function(x,
                   y,
                   itemtypes = c("bernoulli","categorical","gaussian"),
                   penalty = c("lasso","mcp"),
                   nlambda = 100,
                   lambda.max = 3,
                   gamma = 3,
                   lambda = NULL,
                   anchor = NULL,
                   rasch = FALSE,
                   control = list()){

  ##############
  # Preprocess #
  ##############


  #item response warnings
  # if(responsetype == "binary" | any(responsetype == "binary")) {if(any(!(responses == 0 | responses == 1), na.rm = TRUE)) stop("Item responses must be scored 1 for yes/correct and 0 for no/incorrect.", call. = TRUE)}
  if(!(any(itemtypes == "bernoulli") | any(itemtypes == "categorical") | any(itemtypes == "gaussian"))) stop("Item response types must be distributed 'bernoulli', 'categorical', or 'gaussian'.")
  #lambda warnings
  if(any(lambda < 0)) stop("Lambda values must be non-negative.", call. = TRUE)
  if(length(lambda) > 1 & all(diff(lambda) >= 0)) stop("Lambda values must be in descending order (e.g., lambda = c(-2,-1,0)).")
  #anchor warnings
  if(is.null(anchor) & length(lambda) == 1){if(lambda == 0) stop("Anchor item must be specified with lambda = 0.", call. = TRUE)}
  if(!is.null(anchor) & !is.numeric(anchor)) stop("Anchor items must be numeric (e.g., anchor = 1).", call. = TRUE)
  #quadrature warnings
  # if(quadpts < 1) stop("The number of quadrature points must be greater than 1.\n  15 or more points may be required to get accurate estimates.")
  # if(!is.numeric(quadpts)) stop("Quadrature points must be numeric (e.g., quadpts = 15).")

  #data
  responses <- y
  predictors <- x

  #control parameters
  final.control <- list(standardize = TRUE,
                        num_quadpts = 101,
                        tol = 10^-6,
                        maxit = 10000)
  if(length(control) > 0) final.control[names(control)] <- control #user control parameters

  #lambda
  if(is.null(lambda)){
    lambda <- seq(lambda.max**(1/3),0,length.out = nlambda)**3
  }

  #speeds up computation
  responses <- as.matrix(responses)
  predictors <- as.matrix(predictors)
  num_items <- dim(responses)[2]
  samp_size <- dim(responses)[1]
  num_predictors <- dim(predictors)[2]

  #get latent variable values (i.e., predictor values) for quadrature and tracelines
  theta <- seq(-10, 10, length.out = final.control$num_quadpts)

  #turn data into numeric if not already
  if(any(!sapply(responses,function(x)is.numeric(x)))){responses <- sapply(responses,function(x)as.numeric(x))}
  if(any(!sapply(predictors,function(x)is.numeric(x)))){predictors <- sapply(predictors,function(x)as.numeric(x))}

  #get number of itemtypes for number of items
  if(length(itemtypes) == 1){
    itemtypes <- rep(itemtypes, num_items)
  }

  #get item response types
  num_responses <- rep(1,num_items)
  if(any(itemtypes == "bernoulli" | itemtypes == "categorical")){
    responses[,which(itemtypes == "bernoulli" | itemtypes == "categorical")] <-
      apply(responses[,which(itemtypes == "bernoulli" | itemtypes == "categorical")], 2, function(x) as.numeric(as.factor(x)))
    num_responses[which(itemtypes == "bernoulli" | itemtypes == "categorical")] <-
      apply(responses[,which(itemtypes == "bernoulli" | itemtypes == "categorical")], 2, function(x) length(unique(na.omit(x))))
  }

  #standardize predictors
  if(final.control$standardize == TRUE){
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
  final <- list(Lambda = rep(NA,length(lambda)),
                BIC = rep(NA,length(lambda)),
                Impact = replicate(n=length(lambda), NA, simplify = F),
                DIF = replicate(n=length(lambda), NA, simplify = F),
                call = call)

  ##################################
  # Reg-DIF - Loop through lambda #
  ##################################
  message(paste0('Running ',length(lambda),if(length(lambda) == 1){' model.'}else{' models.'}))

  for(pen in 1:length(lambda)){

    #Maximization settings
    lastp <- p
    eps <- Inf
    iter = 0

    ################
    # EM Algorithm #
    ################

    #loop until convergence or maximum number of iterations
    while(eps > final.control$tol & iter < final.control$maxit){

      #E-step: Evaluate Q function with current parameter estimates p
      elist <- Estep_2pl(p,responses,predictors,theta,samp_size,num_items,num_responses,final.control$num_quadpts)

      #M-step: Optimize parameters
      p <- Mstep_2pl_dif(p,responses,predictors,elist,theta,itemtypes,penalty,lambda[pen],gamma,pen,anchor,rasch,final.control$maxit,samp_size,num_responses,num_items,final.control$num_quadpts,num_predictors)

      #Update and check for convergence: Calculate the difference in parameter estimates from current to previous
      eps = sqrt(sum((unlist(p)-unlist(lastp))^2))

      #Update parameter list + B matrix and previous first derivatives
      lastp <- p

      #update the iteration number
      iter = iter + 1
      if(iter == final.control$maxit) warning("EM iteration limit reached without convergence")

      cat('\r',"Models Completed:",pen-1," Iteration:",iter," Change:",round(eps,7)) #print information about optimization
      flush.console()

    } #End EM once converged or reached iteration limit

    ###############
    # Postprocess #
    ###############

    #obtain negative log-likelihood values for all items (** Currently works for only binary items **)
    itemtrace <- lapply(num_responses, function(x) replicate(n=x, matrix(0,nrow=samp_size,ncol=final.control$num_quadpts), simplify = F))
    ll_dif <- ll_dif_pen <- current_pen <- replicate(n=num_items, 0, simplify = F)
    for (item in 1:num_items) { #loop through items
      etable1 <- etable2 <- elist[[1]][[item]]
      etable1[which(!(etable1[,ncol(etable1)] == 1)),] <- 0
      etable2[which(!(etable2[,ncol(etable2)] == 2)),] <- 0
      p_item <- p[[item]]
      etable1 <- etable1[,1:final.control$num_quadpts]
      etable2 <- etable2[,1:final.control$num_quadpts]

      itemtrace[[item]] <- categorical_traceline_pts(p[[item]],theta,predictors,samp_size,num_responses[item],final.control$num_quadpts) #computing probability of endorsement for theta value using current estimate of a and b
      ll_dif[[item]] <- (-1)*(sum(etable2*log(itemtrace[[item]][[2]]), na.rm = TRUE) + sum(etable1*log(itemtrace[[item]][[1]]), na.rm = TRUE))
      current_pen <- (lambda[pen]*sum(c(abs(p_item[grep("c1_itm",names(p_item))]), abs(p_item[grep("a1_itm",names(p_item))])), na.rm = TRUE))
      ll_dif_pen[[item]] <- ll_dif[[item]] - current_pen #Q function we want to minimize
    }
    #obtain likelihood value for latent variable model
    p_impact <- c(p[[num_items+1]],p[[num_items+2]])
    alpha <- predictors %*% p_impact[grep("g",names(p_impact),fixed=T)]
    phi <- exp(predictors %*% p_impact[grep("b",names(p_impact),fixed=T)])
    prior_scores <- t(sapply(1:samp_size, function(x){dnorm(theta, mean = alpha[x], sd = sqrt(phi[x]))}))
    ll_impact <- -1*sum(elist[[2]]*log(prior_scores))

    p2 <- unlist(p)
    ll <- ll_impact + sum(unlist(ll_dif_pen))
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

    final$Lambda[pen] <- lambda[pen]
    final$BIC[pen] <- round(bic,2)
    final$Impact[[pen]] <- round(parms_impact,2)
    final$DIF[[pen]] <- round(parms_dif,2)


    #stop if lambda is too small on first run (this leads to different results b/c of identification constraints on DIF parameters)
    if(is.null(anchor) & pen == 1 & sum(abs(p2[grep(paste0("cov"),names(p2))])) > 0){
      print(coef(final))
      stop("First Lambda value is too small.\n  Two Options:\n  1. Increase first Lambda value large enough to ensure all DIF parameters are removed from the model.\n  2. Provide anchor item(s).", call. = TRUE)
    }

    #stop if there is a large change in DIF parameters
    if(pen > 1){
      second_last <- sum(final$DIF[[pen-1]][,-grep("0",colnames(final$DIF[[pen-1]]))] == 0)
      last <- sum(final$DIF[[pen]][,-grep("0",colnames(final$DIF[[pen-1]]))] == 0)
      if((second_last - last) > (num_predictors*num_items)){
        print(final)
        stop(paste0("Large increase in the number of DIF parameters from iteration ",pen-1," to ",pen,".\n  Two Options:\n  1. Provide smaller differences between lambda values.\n  2. Provide anchor item(s)."), call. = TRUE)
      }
    }

    cat('\r',"Models Completed:",pen," Iteration:",0," Change: -") #print information about optimization
    flush.console()

  } #Terminate Reg-DIF



#Obtain final results
message(paste0('\nFinished.'))
class(final) <- "regDIF"
return(final)

}

