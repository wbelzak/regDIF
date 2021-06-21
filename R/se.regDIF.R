#' Standard Errors for regDIF Model(s)
#'
#' Obtain standard errors for regDIF model(s).
#'
#' @usage
#' se.regDIF(fit,
#'           se.type = "sem",
#'           tau = NULL,
#'           ...)
#'
#' @param fit A regDIF fitted model object of class \code{regDIF}.
#' Upon designating \code{fit}, the default is to obtain standard errors
#' for the best-fitting model according to the minimum BIC model.
#' @param se.type Character value indicating the method of computing standard
#' errors for a regDIF fitted model. Default is "sem", or the supplemental EM
#' algorithm (see Cai, 2008). Other options are in development and not yet
#' supported.
#' @param tau Optional numeric or vector of tau values corresponding to those
#' already fit in \code{fit}.
#' @param ... Additional arguments to pass to regDIF function if different
#' settings are desired.
#'
#' @return Function returns an object of class \code{se.regDIF}
#'
#' @import stats utils
#'
#' @keywords internal
#'
se.regDIF <- function(fit,
                      se.type = "sem",
                      tau = NULL,
                      ...) {

    # Still in development.
    stop("Getting standard errors for regDIF is not yet supported.", call. = FALSE)

    # Obtain data.
    item_data <- fit$data$item.data
    pred_data <- fit$data$pred.data

    # Arguments from regDIF.
    regdif_args <- as.list(fit$call)
    item_type <- if(is.null(regdif_args$item.type)) {
                  NULL
                 } else {
                  regdif_args$item.type
                 }
    pen_type <- if(is.null(regdif_args$pen.type)) {
                  NULL
                } else {
                  regdif_args$pen.type
                }
    tau <- if(is.null(regdif_args$tau)) {
                  fit$tau_vec[[which.min(fit$bic)]]
                } else {
                  regdif_args$tau
                }
    num_tau <- if(is.null(regdif_args$num.tau)) {
                  100
                } else {
                  regdif_args$num.tau
                }
    anchor <- if(is.null(regdif_args$anchor)) {
                  NULL
                } else {
                  regdif_args$anchor
                }
    alpha <- if(is.null(regdif_args$alpha)) {
                  1
                } else {
                  regdif_args$alpha
                }
    gamma <- if(is.null(regdif_args$gamma)) {
                  3
                } else {
                  regdif_args$gamma
                }
    stdz <- if(is.null(regdif_args$stdz)) {
                  TRUE
                } else {
                  regdif_args$stdz
                }
    control <- if(is.null(regdif_args$control)) {
                  list()
                } else {
                  regdif_args$control
                }


    # Pre-process data.
    call <- match.call()
    data_scrub <- preprocess(item_data,
                             pred_data,
                             item_type,
                             pen_type,
                             tau,
                             num_tau,
                             anchor,
                             stdz,
                             control,
                             call)

    # Obtain EM history of minimum BIC model from regDIF model object.
    em_history <- fit$em_history[[which.min(fit$bic)]]

    # Obtain MLEs from EM history.
    mle <- em_history[,ncol(em_history)]

    # Organize MLEs.
    mle_organized <- data_scrub$p
    for(item in 1:data_scrub$num_items) {
      mle_organized[[item]] <- mle[grep(paste0("item",item),names(mle))]
    }
    mle_organized[[data_scrub$num_items+1]] <-
      mle[grep(paste0("mean"),names(mle))]
    mle_organized[[data_scrub$num_items+2]] <-
      mle[grep(paste0("var"),names(mle))]

    # Save MLEs to EM map.
    em_map <- mle_organized

    # Obtain index of SEM "sweet-spot" according to Tian, Cai, et al. (2012).
    obs_ll <- em_history[nrow(em_history),]
    obs_ll_exp_diff <- exp(-diff(obs_ll))
    sem_sweet_spot <- which(obs_ll_exp_diff > .9 & obs_ll_exp_diff < .999)

    # Make space for Jacobian of EM map.
    em_map_jacobian <- last_em_map_jacobian <- eps <- stand_err <- data_scrub$p
    for(item in 1:data_scrub$num_items) {
      em_map_jacobian[[item]][[2]] <-
        last_em_map_jacobian[[item]][[2]] <-
        stand_err[[item]][[2]] <-
        0
    }

    for(parm_vector in 1:(data_scrub$num_items+2)) {
      for(parm in 1:length(eps[[parm_vector]])){
        if(mle_organized[[parm_vector]][[parm]] == 0) next
        eps[[parm_vector]][[parm]] <- Inf
      }
    }

    # Convergence settings.
    iter <- 1
    tol_all <- 1e-2
    tol_parm <- 1e-3
    maxit <- 5000

    while(sum(unlist(eps)) > tol_all && iter < maxit) {

      # Cycle through em iterations in the "sweet-spot".
      for(em_iter in sem_sweet_spot) {

        # Obtain parameter estimates for single em step.
        em_history_current <- em_history[1:(nrow(em_history)-1),em_iter]

        # Cycle through parameter vectors (i.e., item responses and impact
        # equations).
        for(parm_vector in 1:(data_scrub$num_items+2)) {

          # Organize parameter estimates.
          if(parm_vector <= data_scrub$num_items) {
            parm_vector_current <-
              em_history_current[grep(paste0("item",parm_vector),
                                      names(em_history_current))]
          } else if(parm_vector == (data_scrub$num_items+1)){
            parm_vector_current <-
              em_history_current[grep(paste0("mean"),
                                      names(em_history_current))]
          } else {
            parm_vector_current <-
              em_history_current[grep(paste0("var"),
                                      names(em_history_current))]
          }

          # Cycle through parameters within parameter vectors.
          for(parm in 1:length(em_map[[parm_vector]])) {

            # Skip estimates that are 0.
            if(em_map[[parm_vector]][parm] == 0) next
            if(eps[[parm_vector]][parm] < tol_parm) next

            # Update EM map with single parameter from em_history_current.
            em_map[[parm_vector]][[parm]] <- parm_vector_current[[parm]]

            # Single E-step.
            etable <- Estep(em_map,
                            data_scrub$item_data,
                            data_scrub$pred_data,
                            data_scrub$mean_predictors,
                            data_scrub$var_predictors,
                            data_scrub$theta,
                            data_scrub$samp_size,
                            data_scrub$num_items,
                            data_scrub$num_responses,
                            data_scrub$adapt_quad,
                            data_scrub$num_quad)

            # Single M-step.
            if(data_scrub$optim_method == "multi") {
              mout <- Mstep_block(em_map,
                               data_scrub$item_data,
                               data_scrub$pred_data,
                               data_scrub$mean_predictors,
                               data_scrub$var_predictors,
                               etable,
                               data_scrub$item_type,
                               data_scrub$pen_type,
                               data_scrub$tau_vec,
                               alpha,
                               gamma,
                               data_scrub$anchor,
                               data_scrub$samp_size,
                               data_scrub$num_responses,
                               data_scrub$num_items,
                               data_scrub$num_quad,
                               data_scrub$num_predictors)
            } else if(data_scrub$optim_method == "uni") {
              mout <- Mstep_cd(p,
                               data_scrub$item_data,
                               data_scrub$pred_data,
                               data_scrub$mean_predictors,
                               data_scrub$var_predictors,
                               etable,
                               data_scrub$item_type,
                               data_scrub$pen_type,
                               data_scrub$tau_vec,
                               alpha,
                               gamma,
                               data_scrub$anchor,
                               data_scrub$samp_size,
                               data_scrub$num_responses,
                               data_scrub$num_items,
                               data_scrub$num_quad,
                               data_scrub$num_predictors)
            }

            # Get Jacobian of EM map.
            p <- mout$p
            em_map_jacobian_parm <-
              (p[[parm_vector]][[parm]] -
                 mle_organized[[parm_vector]][[parm]]) /
              (parm_vector_current[[parm]] -
                 mle_organized[[parm_vector]][[parm]])

            em_map_jacobian[[parm_vector]][[parm]] <- em_map_jacobian_parm

            # Check for convergence and update.
            eps[[parm_vector]][[parm]] <-
              sqrt((em_map_jacobian[[parm_vector]][[parm]] -
                          last_em_map_jacobian[[parm_vector]][[parm]])**2)
            last_em_map_jacobian <- em_map_jacobian

            # Revert EM map back to MLEs.
            em_map <- mle_organized

            # Update the iteration number.
            iter = iter + 1
            if(iter == maxit) {
              warning("SEM iteration limit reached without convergence")
            }

            # Print SEM iteration
            cat('\r', sprintf("SEM Convergence: Iteration: %d  Change: %f",
                              iter,
                              round(sum(unlist(eps)), nchar(tol_all))))


          }
        }
      }

    }

    # Obtain SEs using EM map.
    complete_info <- fit$complete_ll_info
    for(parm_vector in 1:length(complete_info)) {
      for(parm in 1:length(complete_info[[parm_vector]])) {
        if(em_map_jacobian[[parm_vector]][[parm]] == 0) next
        stand_err[[parm_vector]][[parm]] <-
          sqrt(complete_info[[parm_vector]][[parm]]) *
          (1/(1-em_map_jacobian[[parm_vector]][[parm]]))
      }
    }


}
