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

        names(simtraj) <- names(object@data)[-1]

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

        df
    }
)

setMethod("coef", "fitodeMCMC", function(object) object@coef)

setMethod("vcov", "fitodeMCMC", function(object) object@vcov)

setMethod("stdEr", "fitodeMCMC", function(x) sqrt(diag(vcov(x))))

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

##' show object
##'
##' @param object fitodeMCMC object
##' @docType methods
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
