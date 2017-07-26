##' Plot a fitode object
##' @aliases plot,fitode-method
##' @param x fitode object
##' @param level the confidence level required
##' @param method confidence interval method
##' @param main an overall title for the plot
##' @param xlim the x limit of the plot
##' @param ylim the y limit of the plot
##' @param xlab a label for the x axis
##' @param ylab a label for the y axis
##' @param add (logical) add to an existing plot?
##' @param col.traj colour of the estimated trajectory
##' @param lty.traj line type of the estimated trajectory
##' @param col.conf colour of the confidence intervals
##' @param lty.conf line type of the confidence intervals
##' @param ... additional arguments to be passed on to the plot function
##' @importFrom bbmle plot
##' @docType methods
##' @exportMethod plot
setMethod("plot", signature(x="fitode", y="missing"),
    function(x, level,
             method=c("delta", "mvrnorm", "wmvrnorm"),
             main, xlim, ylim, xlab, ylab, add=FALSE,
             col.traj="black",lty.traj=1,
             col.conf="red",lty.conf=4,
             ...){
        method <- match.arg(method)
        count <- x@data$count
        pred <- predict(x,level,method=method)
        times <- pred[["times"]]
        i.hat <- pred[["mean"]]

        if (missing(main)) main <- "fitode result"
        if (missing(xlab)) xlab <- "time"
        if (missing(ylab)) ylab <- "count"
        if (missing(ylim)) {
            ymin <- min(i.hat, count)
            ymax <- 1.1 * max(i.hat, count)
            ylim <- c(ymin, ymax)
        }
        if (missing(xlim)) xlim <- c(min(times), max(times))

        if (!add) plot(times, count, xlim=xlim, ylim=ylim,
            xlab=xlab, ylab=ylab, main=main, ...)


        lines(times, i.hat, col=col.traj, lty=lty.traj)

        if (!missing(level)) {
            matlines(times, pred[,3:4], col=col.conf, lty=lty.conf)
        }

        invisible()
    }
)


##' Forecast from an ode fit and find confidence interval
##' @param object fitode object
##' @param level the confidence level required
##' @param times time vector to predict over. Default is set to the time frame of the data.
##' @param method confidence interval method. Default is set to Delta method.
##' @param debug print debugging output?
##' @details
##' See vignette for different methods: \code{vignette("details", package="fitode")}
##' @importFrom bbmle predict
##' @importFrom bbmle confint
##' @importFrom MASS mvrnorm
##' @importFrom grDevices adjustcolor
##' @docType methods
##' @exportMethod predict
setMethod("predict", "fitode",
    function(object,
             level,times,
             method=c("delta", "mvrnorm", "wmvrnorm"),
             nsim=1000,
             debug=FALSE){
        if(missing(times)) times <- object@data$times
        method <- match.arg(method)

        model <- object@data$model
        loglik <- object@data$loglik
        parms <- coef(object)

        ss <- solve(model, times, parms, keep_sensitivity=method=="delta")

        expr <- as.expression(object@data$formula[[3]])

        frame <- c(parms, ss@solution)

        mean <- eval(expr, frame)

        df <- data.frame(times=times,mean=mean)

        ## wquant from King et al.
        wquant <- function (x, weights, probs = c(0.025, 0.975)) {
            idx <- order(x)
            x <- x[idx]
            weights <- weights[idx]
            w <- cumsum(weights)/sum(weights)
            rval <- approx(w,x,probs,rule=1)
            rval$y
        }

        if (!missing(level)) {
            nstate <- length(model@state)
            npar <- length(model@par)

            ll <- (1-level)/2

            if (method != "delta") {
                simtraj <- matrix(NA,nrow=length(times),ncol=nsim)
                simpars <- MASS::mvrnorm(nsim,mu=coef(object),
                                   Sigma=vcov(object))

                for (i in 1:nsim) {
                    ss.tmp <- solve(model, times, simpars[i,], keep_sensitivity=method=="delta")
                    frame.tmp <- c(parms, ss.tmp@solution)
                    simtraj[,i] <- eval(expr, frame.tmp)
                }
            }

            cmat <- switch(method,
                delta={
                    dmds <- lapply(model@state, function(s) Deriv(expr, s))
                    dmdp <- lapply(model@par, function(p) Deriv(expr, p))

                    sens <- vector('list', nstate)
                    for(i in 1:nstate) {
                        sens[[i]] <- eval(dmds[[i]], frame) * ss@sensitivity[[i]]
                    }

                    sens_p <- sapply(dmdp, eval, frame)

                    sens <- do.call("+", sens)

                    if(class(sens_p) == "list")
                        sens <- sens + do.call("cbind", sens_p)

                    xvcov <- object@vcov[1:npar,1:npar]
                    if(any(diag(xvcov < 0)))
                        warning("At least one entries in diag(vcov) is negative. Confidence interval may not be accurate.")

                    mvcov <- sens %*% xvcov %*% t(sens)
                    merr <- sqrt(diag(mvcov))
                    z <- -qnorm(ll)
                    cmat <- data.frame(mean - z * merr, mean + z * merr)
                    cmat
                },
                mvrnorm={
                    cmat <- t(apply(simtraj,1,quantile,c(ll,1-ll)))
                    cmat
                },
                wmvrnorm={
                    traj.logLik <- rep(NA, nsim)
                    observation <- object@data$observation

                    for(i in 1:nsim) {
                        traj.logLik[i] <- sum(Eval(loglik, observation, simtraj[,i], simpars[i,]))
                    }

                    ##FIXME: vcov not symmetric for low tolerance?
                    i <- 10
                    while(!isSymmetric(round(vcov(object), i))) i <- i - 1

                    sample.logLik <- mvtnorm::dmvnorm(simpars, coef(object), round(vcov(object), i), log=TRUE)
                    ww <- exp(traj.logLik-sample.logLik)
                    cmat <- t(apply(simtraj, 1, wquant, weights=ww, probs=c(ll, 1-ll)))
                    cmat
                })

            cmat <- setNames(as.data.frame(cmat), c(paste(100*ll, "%"), paste(100*(1-ll), "%")))

            if (debug && method != "delta") {
                matplot(times, simtraj, type="l",col=adjustcolor("black", alpha.f=0.1), lty=1)
                matlines(times, cmat, col=2, lty=1, lwd=2)
            }

            df <- cbind(df, cmat)
        }
        df
    }
)

