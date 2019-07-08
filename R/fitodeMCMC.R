propfun <- function(chol) {
    c(rnorm(ncol(chol), mean=0, sd=1) %*% chol)
}

##' Fit ordinary differential equations model using MCMC
##' @rdname fitodeMCMC
##' @name fitodeMCMC
##' @param model ode model
##' @param data data frame with time column and observation column
##' @param start named vector of starting parameter values
##' @param tcol time column
##' @param vcov
##' @param prior
##' @param chains
##' @param iter
##' @param burnin
##' @param thin
##' @param refresh
##' @param prior.only (logical) sample from prior distribution only?
##' @param link named vector or list of link functions for model parameters
##' @param fixed named vector or list of model parameters to fix
##' @param solver.opts options for ode integration. See \code{\link{ode}}
##' @param solver ode solver
##' @param debug print debugging output?
##' @import coda
##' @export fitodeMCMC
fitodeMCMC <- function(model, data,
                       start, tcol="times",
                       vcov,
                       prior=list(),
                       chains=1, iter=2000, burnin=1000, thin=1,
                       refresh=max(iter/10, 1),
                       prior.only=FALSE,
                       link,
                       fixed=list(),
                       control=list(maxit=1e5),
                       solver.opts=list(method="rk4"),
                       solver=ode,
                       skip.hessian=FALSE,
                       force.hessian=FALSE,
                       use.ginv=TRUE,
                       ...) {
    if (missing(start)) stop("starting parameters must be specified via 'start'")

    if (missing(vcov)) stop("variance covariance matrix of the proposal distribution must be specified via 'vcov'")

    ## TODO: check vcov structure

    if (length(fixed) > 0) model <- fixpar(model, fixed)

    ## turn off sensitivity equations for faster computation
    model <- Transform(model, keep_sensitivity=FALSE)

    modelpar <- model@par

    link <- check_link(model, link)

    if (any(is.na(match(modelpar, names(start))))) {
        stop(
            paste0("`start` must specify the following parameters:\n",
                   "\node parameters: ", paste(model@par, collapse = ", ")
            )
        )
    }

    ## check prior
    if (class(prior)=="list") {
        if (length(prior) == 0) {
            warning("prior distributions must be specified via 'prior'")
            priorlist <- list()
        } else {
            priorlist <- make_prior(model, unlist(link), prior, prior.density=TRUE, model@keep_sensitivity)
        }
    } else if (class(prior)=="function") {
        stop("Not supported yet")
    } else {
        stop("'prior' must be a list of formulas or a function returning log-prior density")
    }

    ## order parameters ...
    start <- start[modelpar]

    link_data <- lapply(link, make.link)

    linklist <- lapply(c("linkfun", "linkinv", "mu.eta"),
                       function(x) lapply(link_data, "[[", x))

    names(linklist) <- c("linkfun", "linkinv", "mu.eta")

    newpar <- Map(function(x, y) ifelse(x=="identity", y, paste(x, y, sep=".")), x=link, y=modelpar)
    newpar <- unname(unlist(newpar))

    names(linklist$linkfun) <- names(linklist$mu.eta) <- newpar

    start <- apply_link(start, linklist, "linkfun")

    dataname <- sapply(lapply(model@observation, "[[", 2), as.character)

    data <- data[,c(tcol, dataname)]

    names(data)[1] <- "times"

    ## returns log-likelihood (instead of negative log-likelihood)
    objfun <- function(model, par, data, solver.opts, solver, linklist, priorlist, prior.only) {
        origpar <- apply_link(par, linklist, "linkinv")
        derivpar <- apply_link(par, linklist, "mu.eta")

        v <- ifelse(prior.only, 0, try(-logLik.sensitivity(origpar, model, data, solver.opts, solver), silent=TRUE))

        if (length(priorlist) > 0) {
            logp <- eval(priorlist$prior.density, as.list(par))
        } else {
            logp <- logpgrad <- 0
        }

        if (inherits(v, "try-error")) {
            return(NA)
        } else {
            ll <- v[1] + logp

            return(ll)
        }
    }

    reslist <- lplist <- vector('list', chains)

    for (nchain in 1:chains) {
        if (refresh > 0)
            message(paste0("MCMC iterations: ", 1, "/", iter, " (Chain ", nchain, ")"))
        mcmcmat <- matrix(NA, nrow=iter, ncol=length(modelpar))
        lpvec <- rep(NA, iter)
        colnames(mcmcmat) <- names(start)

        mcmcmat[1,] <- start
        lpvec[1] <- objfun(model, start, data, solver.opts, solver, linklist, priorlist, prior.only)

        if (iter > 1) {
            cc <- chol(vcov)

            for (i in 2:iter) {
                if (refresh > 0) {
                    if(((i-1) %% refresh)==0) {
                        message(paste0("MCMC iterations: ", i, "/", iter, " (Chain ", nchain, ")"))
                    } else if (i==iter) {
                        message(paste0("MCMC iterations: ", iter, "/", iter, " (Chain ", nchain, ")"))
                    }
                }

                old.theta <- mcmcmat[i-1,]
                new.theta <- old.theta + propfun(cc)
                new.lp <- objfun(model, new.theta, data, solver.opts, solver, linklist, priorlist, prior.only)

                ## this is OK because the proposal distribution is symmetric
                alpha <- exp(new.lp -lpvec[i-1])

                if (!is.finite(alpha)) alpha <- 0

                if (runif(1) < alpha) {
                    mcmcmat[i,] <- new.theta
                    lpvec[i] <- new.lp
                } else {
                    mcmcmat[i,] <- mcmcmat[i-1,]
                    lpvec[i] <- lpvec[i-1]
                }
            }
        }

        ## somewhat based on
        ## https://github.com/LaplacesDemonR/LaplacesDemon/blob/master/R/Thin.R
        keeprows <- burnin + which(rep(1:thin, len=iter-burnin) == 1)

        ## TODO: create a smaller matrix and store only what you need
        ## don't store everything and thin after

        reslist[[nchain]] <- coda::mcmc(mcmcmat[keeprows,], start=burnin+1, end=iter, thin=thin)
        lplist[[nchain]] <- coda::mcmc(lpvec[keeprows], start=burnin+1, end=iter, thin=thin)
    }

    list(
        reslist, lplist
    )
}
