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
        observation <- x@mle2@data$observation
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
        if(missing(times)) times <- object@data$times
        method <- match.arg(method)

        model <- object@model
        loglik <- object@loglik
        parms <- coef(object)

        ## TODO: define a new model without sensitivity if method is not delta

        ss <- ode.solve(model, times, parms,
                        solver.opts=object@mle2@data$solver.opts,
                        solver=object@mle2@data$solver)

        expr <- object@mle2@data$expr

        frame <- c(parms, ss@solution)

        mean <- eval(expr, frame)

        df <- data.frame(times=times[1:length(mean)],mean=mean)

        if (!missing(level)) {
            nstate <- length(model@state)
            npar <- length(model@par)
            linklist <- object@mle2@data$linklist

            ll <- (1-level)/2

            if (method != "delta") {
                vv <- vcov(object, "fitted")
                vv[lower.tri(vv)] <- t(vv)[lower.tri(vv)]

                simtraj <- matrix(NA,nrow=length(mean),ncol=nsim)
                simpars <- MASS::mvrnorm(nsim,mu=coef(object, "fitted"),
                                   Sigma=vv)
                simpars_orig <- t(apply(simpars, 1, apply_link, linklist, "linkinv"))

                for (i in 1:nsim) {
                    ss.tmp <- ode.solve(object@model, times, simpars_orig[i,],
                                        solver.opts=object@mle2@data$solver.opts,
                                        solver=object@mle2@data$solver)
                    frame.tmp <- c(parms, ss.tmp@solution)
                    simtraj[,i] <- eval(expr, frame.tmp)
                }

            }

            cmat <- switch(method,
                delta={
                    sens <- ode.sensitivity(expr,
                                            object@mle2@data$expr.sensitivity,
                                            model,
                                            parms,
                                            object@data$times,
                                            solver.opts=object@mle2@data$solver.opts,
                                            solver=object@mle2@data$solver)$sensitivity

                    fitted_parms <- coef(object, "fitted")

                    mu.eta <- apply_link(fitted_parms, linklist, "mu.eta")[1:npar]

                    sens <- t(t(sens) * mu.eta)

                    fitted.vcov <- vcov(object, "fitted")[1:npar,1:npar]
                    if(any(diag(fitted.vcov < 0)))
                        warning("At least one entries in diag(vcov) is negative. Confidence interval may not be accurate.")

                    mean.vcov <- sens %*% fitted.vcov %*% t(sens)
                    mean.err <- sqrt(diag(mean.vcov))
                    z <- -qnorm(ll)
                    cmat <- data.frame(mean - z * mean.err, mean + z * mean.err)
                    cmat
                },
                mvrnorm={
                    cmat <- t(apply(simtraj,1,quantile,c(ll,1-ll)))
                    cmat
                },
                wmvrnorm={
                    ## wquant from King et al.
                    wquant <- function (x, weights, probs = c(0.025, 0.975)) {
                        idx <- order(x)
                        x <- x[idx]
                        weights <- weights[idx]
                        w <- cumsum(weights)/sum(weights)
                        rval <- approx(w,x,probs,rule=1)
                        rval$y
                    }

                    traj.logLik <- rep(NA, nsim)
                    observation <- object@mle2@data$observation

                    for(i in 1:nsim) {
                        traj.logLik[i] <- sum(Eval(loglik, observation, simtraj[,i], simpars_orig[i,-c(1:npar)]))
                    }

                    sample.logLik <- mvtnorm::dmvnorm(simpars, coef(object, "fitted"), vv, log=TRUE)
                    ww <- exp(traj.logLik-sample.logLik)
                    ww[is.nan(ww) | is.na(ww)] <- 0
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
##' @param scale scale of parameter to be returned
##' @importFrom bbmle coef
##' @docType methods
##' @exportMethod coef
setMethod("coef", "fitode",
    function(object,scale=c("original", "fitted")){
        scale <- match.arg(scale)
        switch(scale,
            original=object@coef,
            fitted=object@mle2@coef
        )
    }
)

##' Extract covariance matrix of a fit
##'
##' @param object fitode object
##' @param scale scale of parameter to be returned
##' @importFrom bbmle vcov
##' @docType methods
##' @exportMethod vcov
setMethod("vcov", "fitode",
    function(object,scale=c("original", "fitted")){
        scale <- match.arg(scale)
        switch(scale,
            original=object@vcov,
            fitted=object@mle2@vcov
        )
    }
)


##' Calculate standard error
##' @importFrom bbmle stdEr
##' @param x fitode object
##' @param scale scale of parameter to be returned
##' @docType methods
##' @exportMethod stdEr
setMethod("stdEr", "fitode", function(x,scale=c("original", "fitted")){sqrt(diag(vcov(x, scale)))})

##' Extract log-likelihood of a fit
##'
##' @param object fitode object
##' @docType methods
##' @exportMethod logLik
setMethod("logLik", "fitode", function(object){-object@min})

##' Profile object
##'
##' @param fitted fitted model object
##' @importFrom bbmle profile
##' @exportMethod profile
setMethod("profile", "fitode",
    function(fitted,
             trace=FALSE) {
        m <- fitted@mle2

        prof <- profile(m, continuation="naive", trace=trace)
        prof
    }
)

##' @importFrom bbmle summary
setMethod("summary","fitode",
    function(object, scale=c("original", "fitted")) {
        scale <- match.arg(scale)
        ss <- summary(object@mle2)

        ## TODO: finish this

    }

)

##' show object
##'
##' @param object fitode object
##' @docType methods
##' @exportMethod show
setMethod("show", "fitode",
    function(object) {
        cat("Model:", object@model@name, "\n")
        cat("\nObservations:\n")
        for(i in 1:length(object@model@observation)) {
            cat(deparse(object@model@observation[[i]]), "\n")
        }

        cat("\nCoefficients:\n")
        print(coef(object))
        cat("\nLog-Likelihood:")
        cat(round(as.numeric(logLik(object)),2),"\n")
        cat("\nlink: ")
        if (length(object@link)==0) {
            cat("none")
        } else {
            cat(names(object@mle2@coef)[!(names(object@mle2@coef) %in% names(object@coef))])
        }
    }
)
