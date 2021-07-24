#' Plot function for regDIF function
#'
#' @param x Fitted regDIF model object.
#' @param y Unused for plotting regDIF model object.
#' @param method Fit statistic to use for identifying DIF effects in plot.
#' @param color.seed Random seed to sample line colors and line types for
#' DIF effects in plot.
#' @param legend.plot Logical indicating whether to plot a legend. Default is \code{TRUE}.
#' @param ... Additional arguments to be passed through to \code{plot}.
#'
#' @rdname plot.regDIF
#'
#' @importFrom graphics abline legend lines plot
#'
#' @return a \code{"plot"} object for a \code{"regDIF"} fit
#'
#' @export
#'
plot.regDIF <-
  function(x, y = NULL, method = "bic", color.seed = 123, legend.plot = TRUE, ...) {

    tau <- log(x$tau_vec)
    if(length(tau) < 2) stop(
      paste0("Must run multiple tau values to plot."),
      call. = TRUE)
    dif.parms <- x$dif[grep(paste0(c("int","slp"),
                                                   collapse = "|"),
                                            rownames(x$dif)), ]
    min.tau <- tau[which.min(unlist(x[method]))]
    dif.min.tau <- dif.parms[,which(tau == min.tau)]
    nonzero.dif <- dif.min.tau[!(dif.min.tau == 0)]

    # # Find first tau with non-zero dif parms.
    #   for(j in 1:ncol(dif.parms)){
    #     if(sum(abs(dif.parms[,j])) > 0){
    #       first.tau <- j
    #       break
    #     } else{
    #       next
    #     }
    #   }

    plot(tau,
         rep(0,length(tau)),
         type = 'l',
         xlim = c(max(tau, na.rm = T), min(tau, na.rm = T)),
         ylim = c(min(dif.parms, na.rm = TRUE), max(dif.parms, na.rm = TRUE)),
         main = "Regularization Path",
         xlab = expression(log(tau)),
         ylab = "Estimate")
    abline(v = min.tau, lty = 2)

      dif.lines <- matrix(NA,ncol=2,nrow=nrow(dif.parms))
      for(i in 1:nrow(dif.parms)){
        if(rownames(dif.parms)[i] %in% names(nonzero.dif)){
          set.seed(color.seed+i)
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
    lines(c(tau, tau[1]),
          rep(0,length(tau)+1),
          type = 'l',
          xlim = c(tau[1],max(tau, na.rm = TRUE)))
    nonzero.dif.lines <- dif.lines[!(dif.min.tau == 0)]
    if(legend.plot) {
      legend("topleft",
             legend = names(nonzero.dif),
             col = nonzero.dif.lines[1:length(nonzero.dif)],
             lty = as.numeric(
               nonzero.dif.lines[(length(nonzero.dif)+1):(length(nonzero.dif)*2)]
             ),
             lwd = 2,
             cex = 0.75,
             bty = "n")
    }


  }
