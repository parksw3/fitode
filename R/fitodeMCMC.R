propfun <- function(chol) {
    c(rnorm(ncol(chol), mean=0, sd=1) %*% chol)
}

##' This function fits ordinary differential equations models to a uni- or
##' multi-variate time series by MCMC using the Metropolis-Hastings update
##' rule. It searches through the parameter space
##' on link scales, which can provide more efficient posterior sampling.
##'
##' @title Fit ordinary differential equations model using MCMC
##' @rdname fitodeMCMC
##' @name fitodeMCMC
##' @param model ode model
##' @param data data frame with time column and observation column
##' @param start named vector of starting parameter values
##' @param tcol time column
##' @param proposal.vcov variance-covariance matrix of a multivariate normal proposal distribution
##' @param prior list of formulas specifying prior distributions
##' @param chains (numeric) number of chains
##' @param iter (numeric) number of iterations per chain
##' @param burnin (numeric) number of burnin interations
##' @param thin (numeric) thining interval between consecutive observations
##' @param refresh (numeric) refresh interval
##' @param prior.only (logical) sample from prior distribution only?
##' @param link named vector or list of link functions for model parameters
##' @param fixed named vector or list of model parameters to fix and their values
##' @param solver.opts options for ode integration. See \code{\link{ode}}
##' @param solver ode solver
##' @param ... additional arguments (unused)
##' @return An object of class ``fitodeMCMC'' as described in \code{\link{fitodeMCMC-class}}.
##' @import coda
##' @importFrom stats runif median cov as.formula
##' @importFrom utils head tail
##' @export fitodeMCMC
fitodeMCMC <- function(model, data,
                       start, tcol="times",
                       proposal.vcov,
                       prior=list(),
                       chains=1, iter=2000, burnin=iter/2, thin=1,
                       refresh=max(iter/10, 1),
                       prior.only=FALSE,
                       link,
                       fixed=list(),
                       solver.opts=list(method="rk4"),
                       solver=ode,
                       ...) {
    call <- match.call()

    if (missing(start))
        stop("starting parameters must be specified via 'start'")

    if (missing(proposal.vcov))
        stop("variance covariance matrix of the proposal distribution must be specified via 'proposal.vcov'")

    if (length(fixed) > 0) model <- fixpar(model, fixed)

    ## turn off sensitivity equations for faster computation
    model <- Transform(model, keep_sensitivity=FALSE)

    link <- check_link(model, link)

    ## check prior
    if (inherits(prior, "list")) {
        if (length(prior) == 0)
 {
            warning("prior distributions must be specified via 'prior'")
            priorlist <- list()
        } else {
            priorlist <- make_prior(model, unlist(link), prior, prior.density=TRUE, keep_grad=FALSE)
        }
    # }  else if (inherits(prior, "function")) {
    #     stop("Not supported yet")
    } else {
        stop("'prior' must be a list of formulas")
    }

    modelpar <- model@par

    ## order parameters ...
    start <- start[modelpar]

    if (any(is.na(match(modelpar, names(start))))) {
        stop(
            paste0("'start' must be a named vector specifying initial values for following parameters:\n",
                   "\node parameters: ", paste(model@par, collapse = ", ")
            )
        )
    }

    link_data <- lapply(link, make.link)

    linklist <- lapply(c("linkfun", "linkinv", "mu.eta"),
                       function(x) lapply(link_data, "[[", x))

    names(linklist) <- c("linkfun", "linkinv", "mu.eta")

    newpar <- Map(function(x, y) ifelse(x=="identity", y, paste(x, y, sep=".")), x=link, y=modelpar)
    newpar <- unname(unlist(newpar))

    names(linklist$linkfun) <- names(linklist$mu.eta) <- newpar

    start <- apply_link(start, linklist, "linkfun")

    names(data)[match(tcol, names(data))] <- "times"

    ## check vcov structure
    ## need to check this after parameters have been transformed
    ## and after fixed parameters have been incorporated
    proposal.msg <- paste0(
        "'proposal.vcov' must be a symmetric (named or unnamed) matrix with following names:\n",
        "\n", paste(newpar, collapse=", ")
    )

    if (!isSymmetric(proposal.vcov)) stop(proposal.msg)

    if (ncol(proposal.vcov) > length(newpar) || nrow(proposal.vcov) > length(newpar)) stop(proposal.msg)

    if (length(colnames(proposal.vcov)) == 0 || length(colnames(proposal.vcov)) == 0) {
        colnames(proposal.vcov) <- rownames(proposal.vcov) <- newpar
    } else {
        if (any(is.na(match(newpar, colnames(proposal.vcov)))) || any(is.na(match(newpar, colnames(proposal.vcov))))) {
            stop(proposal.msg)
        } else {
            proposal.vcov <- proposal.vcov[newpar, newpar]
        }
    }

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

    proposal.chol <- chol(proposal.vcov)

    for (nchain in 1:chains) {
        if (refresh > 0)
            message(paste0("MCMC iterations: ", 1, "/", iter, " (Chain ", nchain, ")"))

        ## somewhat based on
        ## https://github.com/LaplacesDemonR/LaplacesDemon/blob/master/R/Thin.R
        keeprows <- burnin + which(rep(1:thin, len=iter-burnin) == thin)

        mcmcmat <- matrix(NA, nrow=length(keeprows), ncol=length(modelpar))
        lpvec <- rep(NA, length(keeprows))
        colnames(mcmcmat) <- names(start)

        old.theta <- start
        old.lp <- objfun(model, start, data, solver.opts, solver, linklist, priorlist, prior.only)

        if (iter > 1) {
            count <- 1

            for (i in 2:iter) {
                if (refresh > 0) {
                    if(((i-1) %% refresh)==0) {
                        message(paste0("MCMC iterations: ", i, "/", iter, " (Chain ", nchain, ")"))
                    } else if (i==iter) {
                        message(paste0("MCMC iterations: ", iter, "/", iter, " (Chain ", nchain, ")"))
                    }
                }

                new.theta <- old.theta + propfun(proposal.chol)
                new.lp <- objfun(model, new.theta, data, solver.opts, solver, linklist, priorlist, prior.only)

                ## this is OK because the proposal distribution is symmetric
                alpha <- exp(new.lp - old.lp)

                if (!is.finite(alpha)) alpha <- 0

                if (runif(1) < alpha) {
                    old.theta <- new.theta
                    old.lp <- new.lp
                }

                if (i %in% keeprows) {
                    mcmcmat[count,] <- old.theta
                    lpvec[count] <- old.lp
                    count <- count + 1
                }
            }
        }

        ## TODO (maybe): if (scale=="link")
        mcmcmat <- t(apply(mcmcmat, 1, apply_link, linklist, "linkinv"))

        reslist[[nchain]] <- coda::mcmc(mcmcmat, start=head(keeprows,1), end=tail(keeprows,1), thin=thin)
        lplist[[nchain]] <- coda::mcmc(lpvec, start=head(keeprows,1), end=tail(keeprows,1), thin=thin)
    }

    ## return median estimate
    coef <- apply(do.call("rbind", reslist), 2, median)

    ## covariance on the response scale
    vcov <- cov(do.call("rbind", reslist))

    reslist <- coda::as.mcmc.list(reslist)
    lplist <- coda::as.mcmc.list(lplist)

    new("fitodeMCMC",
        call=call,
        model=model, data=data, coef=coef, vcov=vcov,
        mcmc=reslist,
        lp=lplist,
        link=link,
        fixed=as.list(fixed),
        prior=prior,
        details=list(
            chains=chains,
            iter=iter,
            burnin=burnin,
            thin=thin,
            prior.only=prior.only,
            proposal.vcov=proposal.vcov,
            solver.opts=solver.opts,
            solver=solver,
            linklist=linklist
        )
    )
}
