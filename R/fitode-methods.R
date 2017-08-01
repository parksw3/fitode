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
##' @param nsim number of simulations for mvrnorm, wmvrnorm methods
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
             nsim=1000,
             ...){
        method <- match.arg(method)
        observation <- x@data$observation
        pred <- predict(x,level,method=method, nsim=nsim)
        times <- pred[["times"]]
        mean <- pred[["mean"]]

        if (missing(main)) main <- "fitode result"
        if (missing(xlab)) xlab <- "time"
        if (missing(ylab)) ylab <- "observation"
        if (missing(ylim)) {
            ymin <- min(mean, observation)
            ymax <- 1.1 * max(mean, observation)
            ylim <- c(ymin, ymax)
        }
        if (missing(xlim)) xlim <- c(min(times), max(times))

        if (!add) plot(times, observation, xlim=xlim, ylim=ylim,
            xlab=xlab, ylab=ylab, main=main, ...)


        lines(times, mean, col=col.traj, lty=lty.traj)

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
##' @param nsim number of simulations for mvrnorm, wmvrnorm methods
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
             nsim=1000){
        m <- object@mle2

        if(missing(times)) times <- object@data$times
        method <- match.arg(method)

        model <- m@data$model
        loglik <- m@data$loglik
        parms <- coef(m)

        ss <- solve(model, times, parms, keep_sensitivity=method=="delta")

        expr <- as.expression(m@data$formula[[3]])

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
                simpars <- MASS::mvrnorm(nsim,mu=coef(m),
                                   Sigma=vcov(m))

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

                    xvcov <- m@vcov[1:npar,1:npar]
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
                    observation <- m@data$observation

                    for(i in 1:nsim) {
                        traj.logLik[i] <- sum(Eval(loglik, observation, simtraj[,i], simpars[i,-c(1:npar)]))
                    }

                    ##FIXME: vcov not symmetric for low tolerance?
                    i <- 10
                    while(!isSymmetric(round(vcov(m), i))) i <- i - 1

                    sample.logLik <- mvtnorm::dmvnorm(simpars, coef(m), round(vcov(m), i), log=TRUE)
                    ww <- exp(traj.logLik-sample.logLik)
                    cmat <- t(apply(simtraj, 1, wquant, weights=ww, probs=c(ll, 1-ll)))
                    cmat
                })

            cmat <- setNames(as.data.frame(cmat), c(paste(100*ll, "%"), paste(100*(1-ll), "%")))

            df <- cbind(df, cmat)
        }
        df
    }
)

##' Extract parameter of a fit
##' @param object fitode object
##' @param type type of parameter to be returned
##' @importFrom bbmle coef
##' @importFrom bbmle vcov
##' @docType methods
##' @exportMethod coef
setMethod("coef", "fitode",
    function(object,type=c("original", "fitted")){
        type <- match.arg(type)
        cc <- object@mle2@coef
            switch(type,
                original=transpar(cc,
                    object@transforms$transform,
                    object@transforms$inverse) ,
                fitted=cc
            )
        }
)