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
             xlabs, ylabs=nm,
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

            lines(pred.df$times, pred.df$mean, col=col.traj, lty=lty.traj)

            if (!missing(level)) {
                matlines(pred.df$times, pred.df[,3:4], col=col.conf, lty=lty.conf)
            }
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
             method=c("delta", "wmvrnorm"),
             nsim=1000){
        ## TODO: prediction intervals
        if(missing(times)) times <- sort(unique(object@data$times))
        method <- match.arg(method)

        model <- object@model
        parms <- coef(object)
        fixed <- object@fixed

        ## TODO: define a new model without sensitivity if method is not delta

        if (method != "delta" || missing(level)) {
            model <- Transform(model, keep_sensitivity=FALSE)
        }

        ss <- ode.solve(model, times, parms,
                        solver.opts=object@mle2@data$solver.opts,
                        solver=object@mle2@data$solver)

        expr <- object@model@expr

        frame <- c(parms, ss@solution)

        df <- lapply(expr, function(e) {
            m <- eval(e, frame)
            data.frame(
                times=times,
                mean=m
            )
        })

        names(df) <- names(object@data)[-1]

        if (!missing(level)) {
            nstate <- length(model@state)
            npar <- length(model@par)
            linklist <- object@mle2@data$linklist

            ll <- (1-level)/2

            clist <- switch(method,
                delta={
                    sens <- ode.sensitivity(model,
                                            parms,
                                            times,
                                            solver.opts=object@mle2@data$solver.opts,
                                            solver=object@mle2@data$solver)$sensitivity
                    fitted_parms <- coef(object, "fitted")
                    mu.eta <- apply_link(fitted_parms, linklist, "mu.eta")
                    sens <- lapply(sens, function(s) t(t(s) * mu.eta))

                    fitted.vcov <- vcov(object, "fitted")
                    if(any(diag(fitted.vcov < 0)))
                        warning("At least one entries in diag(vcov) is negative. Confidence interval will be accurate.")

                    tmp <- vector('list', length(expr))

                    for (i in 1:length(expr)) {
                        mean.vcov <- sens[[i]] %*% fitted.vcov %*% t(sens[[i]])
                        mean.err <- sqrt(diag(mean.vcov))
                        z <- -qnorm(ll)
                        tmp[[i]] <- data.frame(df[[i]]$mean - z * mean.err, df[[i]]$mean + z * mean.err)
                    }

                    tmp

                },
                wmvrnorm={
                    wmv <- wmvrnorm(object, nsim=nsim)

                    lapply(wmv$simtraj, function(mat) {
                        t(apply(mat, 1, wquant, weights=wmv$weight, probs=c(ll, 1-ll)))
                    })
            })

            clist <- lapply(clist, function(cmat){
                setNames(as.data.frame(cmat), c(paste(100*ll, "%"), paste(100*(1-ll), "%")))
            })

            df <- Map(cbind, df, clist)
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
             which=1:p,
             alpha=0.05,
             trace=FALSE) {
        ## TODO: make this fancier?

        p <- length(fitted@coef)

        prof <- profile(fitted@mle2, which=which, continuation="naive", trace=trace,
                        alpha=alpha)

        prof
    }
)

setMethod("confint", "fitode",
    function (object, parm, level=0.95,
              method=c("delta", "profile", "wmvrnorm"),
              nsim=1000,
              seed) {

        method <- match.arg(method)
        cc <- coef(object)

        ll <- (1-level)/2

        linklist <- object@mle2@data$linklist

        if (missing(parm)) parm <- names(object@coef)

        if (skip.transformation <- all(is.character(parm))) {
            if (!all(parm %in% object@model@par))
                stop("`parm` does not correspond to model parameters.\n",
                     "`parm` must be a vector of model parameters or list of formulas")

            if (method=="profile") {
                prof <- profile(object, which=match(parm, names(object@coef)),
                                alpha=1-level)

                ci0 <- confint(prof, level=level)

                if (length(parm)==1) ci0 <- t(as.matrix(ci0))

                rownames(ci0) <- names(object@mle2@coef)[match(parm, names(object@coef))]

                ci <- apply(ci0, 2, apply_link, linklist, "linkinv")

                if (length(parm)==1) ci <- matrix(c(ci), nrow=1)

                estimate <- matrix(coef(object)[parm], ncol=1)

                res <- cbind(estimate, ci)

                colnames(res) <- c("estimate", paste(100*ll, "%"), paste(100*(1-ll), "%"))
                rownames(res) <- parm

                return(res)
            }

            parm <- lapply(parm, function(x) {
                ee <- as.name(x)
                as.call(list(as.name('~'), ee, ee))
            })

        } else if (is.list(parm)) {
            if (method=="profile")
                stop("profile is not available for non-model parameter confidence intervals")

            ## TODO: don't allow state variables... it gets complicated
        }

        frame <- as.list(cc)

        expr <- lapply(parm, "[[", 3)

        estimate <- try(sapply(expr, eval, frame))

        if (inherits(estimate, "try-error")) {
            stop("Specified formula(s) cannot be evaluated with estimated parameters")
        }

        estimate <- matrix(estimate, ncol=1)

        if (method=="delta") {
            fitted_parms <- coef(object, "fitted")
            fitted_vcov <- vcov(object, "fitted")

            z <- -qnorm(ll)

            if (skip.transformation) {
                est_err <- sqrt(diag(fitted_vcov))

                lwr <- apply_link(fitted_parms - z * est_err, linklist, "linkinv")
                upr <- apply_link(fitted_parms + z * est_err, linklist, "linkinv")

                res <- cbind(estimate, lwr, upr)
            } else {
                expr_sens <- lapply(parm, function(x) Deriv(x[[3]], names(object@coef)))

                mu.eta <- apply_link(fitted_parms, linklist, "mu.eta")

                sens <- t(sapply(expr_sens, function(x) eval(x, frame) * mu.eta))

                est_vcov <- sens %*% fitted_vcov %*% t(sens)

                est_err <- sqrt(diag(est_vcov))

                res <- cbind(estimate, estimate - z * est_err, estimate + z*est_err)
            }
        } else {
            wmv <- wmvrnorm(object, nsim=nsim, seed=seed)

            samp <- matrix(c(apply(wmv$simpars_orig, 1, function(x) sapply(expr, eval, as.list(x)))), ncol=length(parm), byrow=TRUE)

            res <- cbind(estimate, t(apply(samp, 2, wquant, weights=wmv$weight, prob=c(ll, 1-ll))))

        }

        colnames(res) <- c("estimate", paste(100*ll, "%"), paste(100*(1-ll), "%"))
        rownames(res) <- sapply(parm, function(x) as.character(x[[2]]))

        res

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
