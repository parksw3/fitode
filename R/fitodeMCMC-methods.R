##' Computes estimated trajectories and their credible intervals.
##' The estimated trajectories are obtained by taking the median trajectories
##' from the posterior samples.
##'
##' @title Prediction function for fitodeMCMC objects
##' @param object fitodeMCMC object
##' @param level the credible level required
##' @param times time vector to predict over. Default is set to the time frame of the data.
##' @param simplify (logical) simplify output to return estimated trajectories and their
##' credible intervals? If \code{simplify=FALSE}, all posterior trajectories will be returned
##' @return Estimated trajectories and their credible intervals of the fitodeMCMC object
##' @importFrom stats setNames
##' @docType methods
setMethod("predict", "fitodeMCMC",
    function(object,
             level, times,
             simplify=TRUE){
        if(missing(times)) times <- sort(unique(object@data$times))

        ss <- do.call("rbind", object@mcmc)

        model <- object@model

        simtraj <- vector('list', length(model@expr))

        for (j in 1:length(model@expr)) simtraj[[j]] <- matrix(NA, nrow=length(unique(object@data$times)), ncol=nrow(ss))

        for (i in 1:nrow(ss)) {
            ss.tmp <- ode.solve(model, times, ss[i,],
                                solver.opts=object@details$solver.opts,
                                solver=object@details$solver)
            for (k in 1:length(model@expr)) {
                simtraj[[k]][,i] <- eval(model@expr[[k]], c(ss.tmp@solution, ss[i,]))
            }
        }

        names(simtraj) <- sapply(lapply(model@observation, "[[", 2), deparse)

        if (!simplify) return(simtraj)

        df <- lapply(simtraj, function(x) {
            data.frame(
                times=times,
                estimate=apply(x, 1, median)
            )
        })

        if (!missing(level)) {
            ll <- (1-level)/2

            clist <- lapply(simtraj, apply, 1, quantile, probs=c(ll, 1-ll))

            clist <- lapply(clist, function(cmat){
                setNames(as.data.frame(t(cmat)), c(paste(100*ll, "%"), paste(100*(1-ll), "%")))
            })

            df <- Map(cbind, df, clist)
        }

        attr(df, "loglik.ode") <- sapply(object@model@loglik, function(x) x@name)

        df
    }
)

##' Extracts estimated parameters (median of the marginal posterior distributions)
##'
##' @title Extract model coefficients from fitodeMCMC objects
##' @param object fitodeMCMC object
##' @return The estimated median coefficients of the fitodeMCMC object
##' @docType methods
setMethod("coef", "fitodeMCMC", function(object) object@coef)

##' Calculates variance-covariance matrix from posterior samples
##'
##' @title Extract variance-covariance matrix from fitodeMCMC objects
##' @param object fitodeMCMC object
##' @return The variance-covariance matrix of the fitodeMCMC object
##' @docType methods
setMethod("vcov", "fitodeMCMC", function(object) object@vcov)

##' Calculates standard error by taking the square root of the diagonal of the variance-covariance matrix
##' @title Extract standard error from fitodeMCMC objects
##' @param x fitodeMCMC object
##' @return The standard error of the fitodeMCMC object
##' @docType methods
setMethod("stdEr", "fitodeMCMC", function(x) sqrt(diag(vcov(x))))


##' Calculate credible intervals for model parameters and their transformations
##' from posterior samples.
##'
##' @title Calculate credible intervals from fitodeMCMC objects for model parameters and their transformations
##' @param object fitodeMCMC object
##' @param parm character vector specifying model parameters or list of formuals specifying transformations
##' @param level the credible level required
##' @return The credible intervals of the fitodeMCMC object
##' @docType methods
setMethod("confint","fitodeMCMC",
    function(object, parm, level=0.95) {
        ll <- (1-level)/2

        linklist <- object@details$linklist

        ss <- do.call("rbind", object@mcmc)

        if (missing(parm)) parm <- names(object@coef)

        if (skip.transformation <- all(is.character(parm))) {
            if (!all(parm %in% object@model@par))
                stop("`parm` does not correspond to model parameters.\n",
                     "`parm` must be a vector of model parameters or list of formulas")

            parm <- lapply(parm, function(x) {
                ee <- as.name(x)
                as.call(list(as.name('~'), ee, ee))
            })

        }

        expr <- lapply(parm, "[[", 3)

        expr.eval <- matrix(apply(ss, 1, function(x) sapply(expr, eval, as.list(x))), ncol=length(parm), byrow=TRUE)

        res <- cbind(
            apply(expr.eval, 2, median),
            t(apply(expr.eval, 2, quantile, prob=c(ll, 1-ll)))
        )

        colnames(res) <- c("estimate", paste(100*ll, "%"), paste(100*(1-ll), "%"))
        rownames(res) <- sapply(parm, function(x) as.character(x[[2]]))

        res
    }
)

##' Summarize fitodeMCMC object;
##' returns estimate, standard error, credible intervals, effective sample sizes, and gelman-rubin diagnostic
##'
##' @title Summarize fitodeMCMC object
##' @param object fitodeMCMC object
##' @return The summary of the fitodeMCMC object
##' @seealso \code{\link{effectiveSize}} \code{\link{gelman.diag}}
##' @docType methods
setMethod("summary","fitodeMCMC",
    function(object) {
        mm <- matrix(NA, nrow=length(object@coef), ncol=6)

        rownames(mm) <- names(object@coef)
        colnames(mm) <- c("Estimate", "Std. Error", "l-95% CI", "u-95% CI", "Eff. Sample", "Rhat")

        mm[,1] <- object@coef

        if (all(colnames(object@vcov)==names(object@coef)))
            mm[,2] <- stdEr(object)

        mm[,3:4] <- confint(object)[,-1]

        mm[,5] <- round(coda::effectiveSize(object@mcmc), 0)
        mm[,6] <- round(coda::gelman.diag(object@mcmc)[[1]][,1], 2)

        mm
    }
)

##' Show fitodeMCMC object
##'
##' @title Show fitodeMCMC object
##' @param object fitodeMCMC object
##' @return No return value, called for side effects
##' @docType methods
##' @keywords internal
##' @exportMethod show
setMethod("show", "fitodeMCMC",
    function(object) {
        cat("Model:", object@model@name, "\n")
        cat("\nObservations:\n")
        for (i in 1:length(object@model@observation)) {
            cat(deparse(object@model@observation[[i]]), "\n")
        }
        cat("\nPriors:\n")
        if (length(object@prior)==0) {
            cat("\nNone\n")
        } else {
            for (i in 1:length(object@prior)) {
                cat(deparse(object@prior[[i]]), "\n")
            }
        }

        cat("\nCoefficients:\n")
        print(coef(object))
        cat("\nSamples: ")
        cat(paste0(object@details$chains, " chains, each with iter = ", object@details$iter,
                   "; burnin = ", object@details$burnin, "; thin = ", object@details$thin, "\n"))
        cat("\nlink: ")
        if (length(object@link)==0) {
            cat("none")
        } else {
            cat(paste(paste0(names(object@link), " = ", object@link), collapse="; "))
        }
    }
)
