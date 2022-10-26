whichpred <- function(pred, which) {
    tmpattr <- attr(pred, "loglik.ode")[which]
    pred <- pred[which]
    attr(pred, "loglik.ode") <- tmpattr
    pred
}

##' Plot a fitode object
##' @aliases plot,fitode-method
##' @param x fitode object
##' @param level the confidence level required
##' @param which which to plot
##' @param method confidence interval method
##' @param onepage (logical) print all figures on one page?
##' @param xlim x coordinates range
##' @param ylim y coordinates range
##' @param xlabs a label for the x axis
##' @param ylabs a label for the y axis
##' @param col.traj colour of the estimated trajectory
##' @param lty.traj line type of the estimated trajectory
##' @param col.conf colour of the confidence intervals
##' @param lty.conf line type of the confidence intervals
##' @param add add to another plot?
##' @param nsim number of simulations for mvrnorm, wmvrnorm methods
##' @param data (FIXME)
##' @param ... additional arguments to be passed on to the plot function
##' @return No return value, called for side effects
##' @importFrom bbmle plot
##' @docType methods
##' @exportMethod plot
setMethod("plot", signature(x="fitode", y="missing"),
    function(x, level,
             data,
             which,
             method=c("delta", "impsamp", "wmvrnorm"),
             onepage=TRUE,
             xlim, ylim,
             xlabs, ylabs,
             col.traj="black",lty.traj=1,
             col.conf="black",lty.conf=4,
             add=FALSE,
             nsim=1000,
             ...){
        method <- match.arg(method)

        if (method=="impsamp") method <- "wmvrnorm"

        if (missing(data)) data <- x@data

        pred <- predict(x,level,method=method, nsim=nsim)

        if (missing(which)) which <- 1:length(pred)
        if (is.character(which)) which <- match(which,names(pred))

        pred <- whichpred(pred, which)

        plot_internal(pred, data, onepage, xlim, ylim, xlabs, ylabs, col.traj, lty.traj, col.conf, lty.conf, add, ...)
    }
)

##' Plot a fitodeMCMC object
##' @aliases plot,fitodeMCMC-method
##' @param x fitodeMCMC object
##' @param level the confidence level required
##' @param which which to plot
##' @param onepage (logical) print all figures on one page?
##' @param xlim x coordinates range
##' @param ylim y coordinates range
##' @param xlabs a label for the x axis
##' @param ylabs a label for the y axis
##' @param col.traj colour of the estimated trajectory
##' @param lty.traj line type of the estimated trajectory
##' @param col.conf colour of the confidence intervals
##' @param lty.conf line type of the confidence intervals
##' @param add add to another plot?
##' @param data (FIXME)
##' @param ... additional arguments to be passed on to the plot function
##' @return No return value, called for side effects
setMethod("plot", signature(x="fitodeMCMC", y="missing"),
    function(x, level,
             data,
             which,
             onepage=TRUE,
             xlim, ylim,
             xlabs, ylabs,
             col.traj="black",lty.traj=1,
             col.conf="black",lty.conf=4,
             add=FALSE,
             ...){
        if (missing(data)) data <- x@data

        pred <- predict(x,level,simplify=TRUE)

        if (missing(which)) which <- 1:length(pred)
        if (is.character(which)) which <- match(which,names(pred))

        pred <- whichpred(pred, which)

        plot_internal(pred, data, onepage, xlim, ylim, xlabs, ylabs, col.traj, lty.traj, col.conf, lty.conf, add, ...)
    }
)

##' Internal function for plotting methods
##' @param pred prediction objects
##' @param data observed data
##' @param onepage (logical) print all figures on one page?
##' @param xlim x coordinates range
##' @param ylim y coordinates range
##' @param xlabs a label for the x axis
##' @param ylabs a label for the y axis
##' @param col.traj colour of the estimated trajectory
##' @param lty.traj line type of the estimated trajectory
##' @param col.conf colour of the confidence intervals
##' @param lty.conf line type of the confidence intervals
##' @param add add to another plot?
##' @param ... additional arguments to be passed on to the plot function
##' @importFrom graphics par lines matlines
plot_internal <- function(pred,
                          data,
                          onepage=TRUE,
                          xlim, ylim,
                          xlabs, ylabs,
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

    nm <- names(pred)

    if (missing(ylabs)) ylabs <- nm

    if (missing(xlabs)) xlabs <- rep("times", length(pred))

    auto_ylim <- missing(ylim)
    auto_xlim <- missing(xlim)

    for (i in 1:length(pred)) {
        obs.df <- data.frame(
            x=data$times,
            y=eval(parse(text=nm[i]), data)
        )

        pred.df <- pred[[i]]

        if (attr(pred, "loglik.ode")[i]=="dlnorm")
            pred.df[,-1] <- exp(pred.df[,-1])

        if (auto_ylim) {
            ymin <- 0.95 * min(min(unlist(pred.df[,-1])), obs.df$y, na.rm=TRUE)
            ymax <- 1.05 * max(max(unlist(pred.df[,-1])), obs.df$y, na.rm=TRUE)
            ylim <- c(ymin, ymax)
        }

        if (auto_xlim) {
            xlim <- c(min(obs.df$x), max(obs.df$x))
        }

        if (!add) plot(obs.df, xlim=xlim, ylim=ylim,
                       xlab=xlabs[i], ylab=ylabs[i], ...)

        lines(pred.df$times, pred.df$estimate, col=col.traj, lty=lty.traj)

        if (ncol(pred.df) > 2) {
            matlines(pred.df$times, pred.df[,3:4], col=col.conf, lty=lty.conf)
        }
    }

    invisible()
}
