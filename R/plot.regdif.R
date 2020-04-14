#' Plot function for regDIF function
#'
#' @param x Fitted regDIF model object.
#' @param y Unused for plotting regDIF model object.
#' @param method Fit statistic to use for identifying DIF effects in plot.
#' @param legend.seed Random seed to sample line colors and line types for DIF effects in plot.
#' @param ... Additional arguments to be passed through to \code{plot}.
#'
#' @rdname plot.regDIF
#' @importFrom graphics abline legend lines plot
#' @export
#'

plot.regDIF <-
  function(x, y = NULL, method = "bic", legend.seed = 123, ...) {

    lambda <- x$lambda
    if(length(lambda) < 2) stop(paste0("Must run multiple lambda values to plot."), call. = TRUE)
    dif.parms <- x$dif.item.parms[grep(paste0(c("int","slp"),collapse = "|"),rownames(x$dif.item.parms)),]
    min.lambda <- lambda[which.min(unlist(x[method]))]
    dif.min.lambda <- dif.parms[,which(lambda == min.lambda)]
    nonzero.dif <- dif.min.lambda[!(dif.min.lambda == 0)]
    if(length(nonzero.dif) == 0) stop(paste0("No DIF effects in final model to plot."), call. = TRUE)
    #find first lambda with non-zero dif parms
      for(j in 1:ncol(dif.parms)){
        if(sum(abs(dif.parms[,j])) > 0){
          first.lambda <- j
          break
        } else{
          next
        }
      }

    plot(lambda, rep(0,length(lambda)), type = 'l', xlim = c(lambda[first.lambda]+.5,min(lambda)), main = "Regularization Paths", xlab = expression(lambda), ylab = "Estimate")
    abline(v = min.lambda, lty = 2)

      dif.lines <- matrix(NA,ncol=2,nrow=nrow(dif.parms))
      for(i in 1:nrow(dif.parms)){
        if(rownames(dif.parms)[i] %in% names(nonzero.dif)){
          set.seed(legend.seed+i)
          linecolor <- sample(grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)], 1)
          linetype <- sample(1:6,1)
          lwdnum <- 2
        } else {
          linecolor <- 'gray72'
          linetype <- 1
          lwdnum <- 1
        }

        lines(lambda,dif.parms[i,], col = linecolor, lty = linetype, lwd = lwdnum)
        dif.lines[i,1] <- linecolor
        dif.lines[i,2] <- linetype

      }
    lines(c(lambda[first.lambda]+.5,lambda), rep(0,length(lambda)+1), type = 'l', xlim = c(lambda[first.lambda]+.5,min(lambda)))
    nonzero.dif.lines <- dif.lines[!(dif.min.lambda == 0)]
    legend("topleft", legend = names(nonzero.dif), col = nonzero.dif.lines[1:length(nonzero.dif)], lty = as.numeric(nonzero.dif.lines[(length(nonzero.dif)+1):(length(nonzero.dif)*2)]), lwd = 2, cex = 0.75)

  }
