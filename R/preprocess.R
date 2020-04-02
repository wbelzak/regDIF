##############
# Preprocess #
##############

preprocess <-
  function(x,
           y,
           itemtypes,
           penalty,
           nlambda,
           lambda.max,
           lambda,
           anchor,
           rasch,
           control,
           call){

  #preprocess warnings
  if(!(any(itemtypes == "bernoulli") | any(itemtypes == "categorical") | any(itemtypes == "gaussian"))) stop("Item response types must be distributed 'bernoulli', 'categorical', or 'gaussian'.")
  if(any(lambda < 0)) stop("Lambda values must be non-negative.", call. = TRUE)
  if(length(lambda) > 1 & all(diff(lambda) >= 0)) stop("Lambda values must be in descending order (e.g., lambda = c(-2,-1,0)).")
  if(is.null(anchor) & length(lambda) == 1){if(lambda == 0) stop("Anchor item must be specified with lambda = 0.", call. = TRUE)}
  if(!is.null(anchor) & !is.numeric(anchor)) stop("Anchor items must be numeric (e.g., anchor = 1).", call. = TRUE)

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
  final <- list(Lambda = rep(NA,length(lambda)),
                AIC = rep(NA,length(lambda)),
                BIC = rep(NA,length(lambda)),
                Impact = replicate(n=length(lambda), NA, simplify = F),
                DIF = replicate(n=length(lambda), NA, simplify = F),
                call = call)

  return(list(p = p,
              final = final,
              responses = responses,
              predictors = predictors,
              itemtypes = itemtypes,
              final.control = final.control,
              lambda = lambda,
              theta = theta,
              num_responses = num_responses,
              num_predictors = num_predictors,
              samp_size = samp_size,
              num_items = num_items))

}
