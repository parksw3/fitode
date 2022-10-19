##' Computes estimated trajectories and their confidence intervals (using either
##' the delta method or importance sampling).
##'
##' @title Prediction function for fitode objects
##' @param object fitode object
##' @param level the confidence level required
##' @param times time vector to predict over. Default is set to the time frame of the data.
##' @param method confidence interval method. Default is set to Delta method.
##' @param nsim number of simulations for mvrnorm, wmvrnorm methods
##' @return The estimated trajectories and their confidence intervals of the fitode object
##' @importFrom bbmle predict
##' @importFrom bbmle confint
##' @importFrom MASS mvrnorm
##' @importFrom grDevices adjustcolor
##' @importFrom stats qnorm quantile
##' @docType methods
##' @exportMethod predict
setMethod("predict", "fitode",
    function(object,
             level,times,
             method=c("delta", "impsamp", "wmvrnorm"),
             nsim=1000){
        if (missing(times)) times <- sort(unique(object@data$times))
        method <- match.arg(method)

        if (method=="impsamp") method <- "wmvrnorm"

        model <- object@model
        parms <- coef(object)
        fixed <- object@fixed

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
                estimate=m
            )
        })

        names(df) <- sapply(lapply(model@observation, "[[", 2), deparse)

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
                    fitted_parms <- coef(object, "links")
                    mu.eta <- apply_link(fitted_parms, linklist, "mu.eta")
                    sens <- lapply(sens, function(s) t(t(s) * mu.eta))

                    fitted.vcov <- vcov(object, "links")
                    if(any(diag(fitted.vcov < 0)))
                        warning("At least one entries in diag(vcov) is negative. Confidence interval will be accurate.")

                    tmp <- vector('list', length(expr))

                    for (i in 1:length(expr)) {
                        mean.vcov <- sens[[i]] %*% fitted.vcov %*% t(sens[[i]])
                        mean.err <- sqrt(diag(mean.vcov))
                        z <- -qnorm(ll)
                        tmp[[i]] <- data.frame(df[[i]]$estimate - z * mean.err, df[[i]]$estimate + z * mean.err)
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

        attr(df, "loglik.ode") <- sapply(object@model@loglik, function(x) x@name)

        df
    }
)

##' Extracts estimated parameters (either on response scales or link scales)
##'
##' @title Extract model coefficients from fitode objects
##' @param object fitode object
##' @param type type of coefficients. The default (\code{type=response}) is on the
##' response scale; this is the scale on which the model parameters are defined.
##' Alternatively, \code{type=link} can be used to obtain parameters on the estimated scale.
##' @return The estimated coefficients of the fitode object
##' @importFrom bbmle coef
##' @docType methods
##' @exportMethod coef
setMethod("coef", "fitode",
    function(object,type=c("response", "links")){
        type <- match.arg(type)
        switch(type,
            response=object@coef,
            links=object@mle2@coef
        )
    }
)

##' Extracts variance-covariance matrix (either on response scales or link scales)
##'
##' @title Extract variance-covariance matrix from fitode objects
##' @param object fitode object
##' @param type type of covariance matrix. The default (\code{type=response}) is on the
##' response scale; this is the scale on which the model parameters are defined.
##' Alternatively, \code{type=link} can be used to obtain the covariance matrix on the estimated scale.
##' @return The variance-covariance matrix of the fitode object
##' @importFrom bbmle vcov
##' @docType methods
##' @exportMethod vcov
setMethod("vcov", "fitode",
    function(object,type=c("response", "links")){
        type <- match.arg(type)
        switch(type,
            response=object@vcov,
            links=object@mle2@vcov
        )
    }
)


##' Calculates standard error by taking the square root of the diagonal matrix
##' @title Extract standard error from fitode objects
##' @importFrom bbmle stdEr
##' @param x fitode object
##' @param type type of standard error. The default (\code{type=response}) is on the
##' response scale; this is the scale on which the model parameters are defined.
##' Alternatively, \code{type=link} can be used to obtain standard errors on the estimated scale.
##' @return The standard error of the fitode object
##' @docType methods
##' @exportMethod stdEr
setMethod("stdEr", "fitode", function(x,type=c("response", "links")){sqrt(diag(vcov(x, type)))})

##' Extract log-likelihood of a fit
##' @title Extract log-likelihood
##' @param object fitode object
##' @return The log-likelihood of the fitode object
##' @docType methods
##' @exportMethod logLik
setMethod("logLik", "fitode", function(object){-object@min})

##' Profile fitode objects
##'
##' @param fitted fitted model object
##' @importFrom bbmle profile
##' @param which which parameter(s) to profile? (integer value)
##' @param alpha critical level
##' @param trace trace progress of computations?
##' @param ... additional arguments passed to mle2 profiling method
##' @return The log-likelihood profile of the fitode object
##' @exportMethod profile
setMethod("profile", "fitode",
    function(fitted,
             which=1:p,
             alpha=0.05,
             trace=FALSE,
             ...) {
        ## TODO: make this fancier?
        p <- length(fitted@coef)
        prof <- profile(fitted@mle2, which=which, continuation="naive", trace=trace,
                        alpha=alpha,...)

        prof
    }
)

##' Calculate confidence intervals for model parameters and their transformations using
##' (1) delta method, (2) profile likelihood, and (3) importance sampling.
##'
##' @title Calculate confidence intervals from fitode objects for model parameters and their transformations
##' @param object fitode object
##' @param parm character vector specifying model parameters or list of formuals specifying transformations
##' @param level the confidence level required
##' @param method method for calculating confidence intervals
##' @param nsim number of simulations to be used for importance sampling
##' @param seed seed
##' @param ... extra arguments passed to profiling method
##' @return The confidence intervals for model parameters and their transformations of the fitode object
##' @docType methods
##' @exportMethod confint
setMethod("confint", "fitode",
    function (object, parm, level=0.95,
              method=c("delta", "profile", "impsamp", "wmvrnorm"),
              nsim=1000,
              seed,
              ...) {

        method <- match.arg(method)

        if (method=="impsamp") method <- "wmvrnorm"

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
                                alpha=1-level, ...)

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
                stop("profile is only available for model parameters")

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
            fitted_parms <- coef(object, "links")
            fitted_vcov <- vcov(object, "links")

            z <- -qnorm(ll)

            if (skip.transformation) {
                fitted_parms <- fitted_parms[match(sapply(expr, as.character), object@model@par)]
                est_err <- sqrt(diag(fitted_vcov))[match(sapply(expr, as.character), object@model@par)]

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

##' Summarize fitode objects;
##' returns estimate, standard error, and confidence intervals
##'
##' @title Summarize fitode object
##' @param object fitode object
##' @return The summary of the fitode object
##' @importFrom bbmle summary
##' @docType methods
##' @exportMethod summary
setMethod("summary","fitode",
    function(object) {
        mm <- matrix(NA, nrow=length(object@coef), ncol=4)

        rownames(mm) <- names(object@coef)
        colnames(mm) <- c("Estimate", "Std. Error", "l-95% CI", "u-95% CI")

        mm[,1] <- object@coef

        if (all(colnames(object@vcov)==names(object@coef)))
            mm[,2] <- stdEr(object)

        mm[,3:4] <- confint(object)[,-1]

        mm
    }
)

##' Show fitode objects
##'
##' @title Show fitode objects
##' @param object fitode object
##' @return No return value, called for side effects
##' @docType methods
##' @exportMethod show
##' @keywords internal
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
            cat(paste(paste0(names(object@link), " = ", object@link), collapse="; "))
        }
        cat("\n")
    }
)
