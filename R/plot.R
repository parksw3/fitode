##' Plot a fitode object
##' @aliases plot,fitode-method
##' @param x fitode object
##' @param level the confidence level required
##' @param which which to plot
##' @param method confidence interval method
##' @param onepage (logical) print all figures on one page?
##' @param xlabs a label for the x axis
##' @param ylabs a label for the y axis
##' @param col.traj colour of the estimated trajectory
##' @param lty.traj line type of the estimated trajectory
##' @param col.conf colour of the confidence intervals
##' @param lty.conf line type of the confidence intervals
##' @param add add to another plot?
##' @param nsim number of simulations for mvrnorm, wmvrnorm methods
##' @param ... additional arguments to be passed on to the plot function
##' @importFrom bbmle plot
##' @docType methods
##' @exportMethod plot
setMethod("plot", signature(x="fitode", y="missing"),
    function(x, level,
             which,
             method=c("delta", "wmvrnorm"),
             onepage=TRUE,
             xlabs, ylabs,
             col.traj="black",lty.traj=1,
             col.conf="black",lty.conf=4,
             add=FALSE,
             nsim=1000,
             ...){
        method <- match.arg(method)
        data <- x@data

        pred <- predict(x,level,method=method, nsim=nsim)

        if (missing(which)) which <- 1:length(pred)
        if (is.character(which)) which <- match(which,names(pred))

        pred <- pred[which]

        plot_internal(pred, data, onepage, xlabs, ylabs, col.traj, lty.traj, col.conf, lty.conf, add, ...)
    }
)

setMethod("plot", signature(x="fitodeMCMC", y="missing"),
    function(x, level,
             which,
             onepage=TRUE,
             xlabs, ylabs,
             col.traj="black",lty.traj=1,
             col.conf="black",lty.conf=4,
             add=FALSE,
             ...){
        method <- match.arg(method)
        data <- x@data

        pred <- predict(x,level,simplify=TRUE)

        if (missing(which)) which <- 1:length(pred)
        if (is.character(which)) which <- match(which,names(pred))

        pred <- pred[which]

        plot_internal(pred, data, onepage, xlabs, ylabs, col.traj, lty.traj, col.conf, lty.conf, add, ...)
    }
)

plot_internal <- function(pred,
                          data,
                          onepage=TRUE,
                          xlabs, ylabs=nm,
                          col.traj="black",lty.traj=1,
                          col.conf="black",lty.conf=4,
                          add=FALSE,
                          ...) {
    ## from bbmle
    if (onepage && !add) {
        nplots <- length(pred)
        ## Q: should we reset par(mfrow), or par(mfg), anyway?
        if (prod(par("mfcol")) < nplots) {
            rows <- ceiling(round(sqrt(nplots)))
            columns <- ceiling(nplots/rows)
            mfrow_orig <- par(mfrow=c(rows,columns))

            on.exit(par(mfrow_orig))
        }
    }


    if (missing(ylabs)) ylabs <- names(pred)

    if (missing(xlabs)) xlabs <- rep("times", length(pred))

    for (i in 1:length(which)) {
        obs.df <- data.frame(
            x=data$times,
            y=data[[nm[i]]]
        )

        pred.df <- pred[[i]]

        ymin <- 0.95 * min(min(unlist(pred.df[,-1])), obs.df$y, na.rm=TRUE)
        ymax <- 1.05 * max(max(unlist(pred.df[,-1])), obs.df$y, na.rm=TRUE)
        ylim <- c(ymin, ymax)
        xlim <- c(min(obs.df$x), max(obs.df$x))

        if (!add) plot(obs.df, xlim=xlim, ylim=ylim,
                       xlab=xlabs[i], ylab=ylabs[i], ...)

        lines(pred.df$times, pred.df$estimate, col=col.traj, lty=lty.traj)

        if (ncol(pred.df) > 2) {
            matlines(pred.df$times, pred.df[,3:4], col=col.conf, lty=lty.conf)
        }
    }

    invisible()

}
