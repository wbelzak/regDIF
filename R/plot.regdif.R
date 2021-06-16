#' Plot function for regDIF function
#'
#' @param x Fitted regDIF model object.
#' @param y Unused for plotting regDIF model object.
#' @param method Fit statistic to use for identifying DIF effects in plot.
#' @param legend.seed Random seed to sample line colors and line types for
#' DIF effects in plot.
#' @param ... Additional arguments to be passed through to \code{plot}.
#'
#' @rdname plot.regDIF
#'
#' @importFrom graphics abline legend lines plot
#'
#' @export
#'
plot.regDIF <-
  function(x, y = NULL, method = "bic", legend.seed = 123, ...) {

    tau <- x$tau_vec
    if(length(tau) < 2) stop(
      paste0("Must run multiple tau values to plot."),
      call. = TRUE)
    dif.parms <- x$dif[grep(paste0(c("int","slp"),
                                                   collapse = "|"),
                                            rownames(x$dif)), ]
    min.tau <- tau[which.min(unlist(x[method]))]
    dif.min.tau <- dif.parms[,which(tau == min.tau)]
    nonzero.dif <- dif.min.tau[!(dif.min.tau == 0)]
    if(length(nonzero.dif) == 0) {
      stop(paste0("No DIF effects in final model to plot."), call. = TRUE)
    }
    # Find first tau with non-zero dif parms.
      for(j in 1:ncol(dif.parms)){
        if(sum(abs(dif.parms[,j])) > 0){
          first.tau <- j
          break
        } else{
          next
        }
      }

    plot(tau,
         rep(0,length(tau)),
         type = 'l',
         xlim = c(tau[first.tau]+.1,min(tau)),
         main = "Regularization Paths",
         xlab = expression(tau),
         ylab = "Estimate")
    abline(v = min.tau, lty = 2)

      dif.lines <- matrix(NA,ncol=2,nrow=nrow(dif.parms))
      for(i in 1:nrow(dif.parms)){
        if(rownames(dif.parms)[i] %in% names(nonzero.dif)){
          set.seed(legend.seed+i)
          linecolor <- sample(grDevices::colors()[grep('gr(a|e)y',
                                                       grDevices::colors(),
                                                       invert = T)], 1)
          linetype <- sample(1:6,1)
          lwdnum <- 2
        } else {
          linecolor <- 'gray72'
          linetype <- 1
          lwdnum <- 1
        }

        lines(tau,dif.parms[i,], col = linecolor, lty = linetype, lwd = lwdnum)
        dif.lines[i,1] <- linecolor
        dif.lines[i,2] <- linetype

      }
    lines(c(tau[first.tau]+.1,tau),
          rep(0,length(tau)+1),
          type = 'l',
          xlim = c(tau[first.tau]+.1,min(tau)))
    nonzero.dif.lines <- dif.lines[!(dif.min.tau == 0)]
    legend("topleft",
           legend = names(nonzero.dif),
           col = nonzero.dif.lines[1:length(nonzero.dif)],
           lty = as.numeric(
             nonzero.dif.lines[(length(nonzero.dif)+1):(length(nonzero.dif)*2)]
             ),
           lwd = 2,
           cex = 0.75)

  }
