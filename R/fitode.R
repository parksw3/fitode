##' Set up link functions for model parameters -- for \code{\link{fitode}}
##' internal usage; assumes log-link by default if link functions are
##' not specified
##'
##' @title Set up link functions for model parameters
##' @param link named list or vector of strings specifying link functions
##' @param par (character) model parameters
##' @seealso \code{\link{make.link}}
##' @keywords internal
##' @return list of strings specifying link functions for each model parameter
set_link <- function(link, par) {
    link_default <- as.list(rep("log", length(par)))
    names(link_default) <- par

    if (!missing(link)) link_default[names(link)] <- link

    link_default
}

##' Apply link functions to model parameters -- for \code{\link{fitode}}
##' internal usage. \code{type=linkfun} applies link functions;
##' \code{type=linkinv} applies inverse link functions; \code{type=mu.eta}
##' returns derivative of inverse link functions
##'
##' @title Apply link functions to model parameters
##' @param par named vector of parameter values
##' @param linklist list containing \code{linkfun}, \code{linkinv}, and \code{mu.eta} for each link
##' @param type string specifying which function should be applied
##' @keywords internal
##' @seealso \code{\link{make.link}}
apply_link <- function(par, linklist, type=c("linkfun", "linkinv", "mu.eta")) {
    type <- match.arg(type)

    if (type=="linkinv" || type=="mu.eta") {
        filter <- match(names(par), names(linklist$linkfun))
    } else {
        filter <- match(names(par), names(linklist$linkinv))
    }

    ff <- linklist[[type]][filter]
    pp <- unlist(Map(function(x, fun) fun(x), x=par, fun=ff))
    names(pp) <- names(ff)

    pp
}

##' Check link functions for model parameters -- for \code{\link{fitode}}
##' internal usage
##'
##' @title Check link functions
##' @param model odemodel object
##' @param link named vector specifying link functions
##' @keywords internal
check_link <- function(model, link) {
    if (missing(link)) return(as.list(model@link))

    if (any(is.na(match(names(link), model@par)))) stop("Some link functions do not correspond to the model parameters.")

    link <- set_link(link, model@par)

    link
}

##' Fix parameters of an ODE model to a constant value by transforming
##' the model (using Transform() function)
##'
##' @title Fix parameters of an ODE model
##' @param model odemodel object
##' @param fixed named vector or list of model parameters to fix
##' @keywords internal
fixpar <- function(model, fixed) {
    fixed <- as.list(fixed)
    if (any(!(names(fixed) %in% model@par)))
        stop("`fixed must be a named vector/list whose names correspond to model parameters")

    tlist <- vector('list', length(fixed))

    for (i in 1:length(fixed)) {
        tlist[[i]] <- as.formula(as.call(c(as.symbol("~"), as.symbol(names(fixed)[i]), unname(fixed[i]))))
    }

    par <- model@par[!(model@par %in% names(fixed))]

    model <- Transform(
        model,
        tlist,
        par
    )

    model
}

##' This function fits ordinary differential equations models to a uni- or
##' multi-variate time series by maximum likelihood.
##' It relies on sensitivity equations to compute
##' gradients of the estimated trajectory with respect to model parameters.
##' This allows one to use gradient-based optimization algorithms, which can
##' provide more robust estimation.
##'
##' @title Fit ordinary differential equations model
##' @rdname fitode
##' @name fitode
##' @param model odemodel object
##' @param data data frame with a time column and observation columns
##' @param start named vector of starting parameter values
##' @param tcol (character) time column
##' @param method optimization method
##' @param optimizer optimizer
##' @param link named vector or list of link functions for model parameters
##' @param fixed named vector or list of model parameters to fix and their values
##' @param prior list of formulas specifying prior distributions
##' @param prior.density (logical) should priors represent probability distributions?
##' @param control see \code{\link{optim}}
##' @param solver.opts options for ode integration. See \code{\link{ode}}
##' @param solver ode solver
##' @param skip.hessian skip hessian calculation
##' @param force.hessian (logical) calculate the hessian numerically instead of taking the jacobian of the gradients based on sensitivity equations
##' @param use.ginv (logical) use generalized inverse (\code{\link{ginv}}) to compute approximate vcov
##' @param ... mle2 arguments
##' @import bbmle
##' @importFrom numDeriv jacobian hessian
##' @importFrom MASS ginv
##' @seealso \code{\link{mle2}}
##' @export fitode
fitode <- function(model, data,
                   start, tcol="times",
                   method="BFGS",
                   optimizer="optim",
                   link,
                   fixed=list(),
                   prior=list(),
                   prior.density=TRUE,
                   control=list(maxit=1e5),
                   solver.opts=list(method="rk4"),
                   solver=ode,
                   skip.hessian=FALSE,
                   force.hessian=FALSE,
                   use.ginv=TRUE,
                   ...) {
    if (missing(start)) stop("starting parameters must be specified via `start'")

    if (length(fixed) > 0) model <- fixpar(model, fixed)

    modelpar <- model@par

    link <- check_link(model, link)

    if (any(is.na(match(modelpar, names(start))))) {

        stop(
            paste0("'start' must be a named vector specifying initial values for following parameters:\n",
                "\node parameters: ", paste(model@par, collapse = ", ")
            )
        )
    }

    if (length(prior) == 0) {
        priorlist <- list()
    } else {
        priorlist <- make_prior(model, unlist(link), prior, prior.density, model@keep_sensitivity)
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

    keep_sensitivity <- model@keep_sensitivity

    dataname <- sapply(lapply(model@observation, "[[", 2), as.character)

    data <- data[,c(tcol, dataname)]

    names(data)[1] <- "times"

    dataarg <- list(model=model, data=data, solver.opts=solver.opts, solver=solver, linklist=linklist,
                    priorlist=priorlist)

    f.env <- new.env()
    ## set initial values
    assign("oldnll",NULL,f.env)
    assign("oldpar",NULL,f.env)
    assign("oldgrad",NULL,f.env)

    objfun <- function(par, data, solver.opts, solver, linklist, priorlist) {
        if (identical(par,oldpar)) {
            return(oldnll)
        }
        origpar <- apply_link(par, linklist, "linkinv")
        derivpar <- apply_link(par, linklist, "mu.eta")

        v <- try(logLik.sensitivity(origpar, model, data, solver.opts, solver), silent=TRUE)

        if (length(priorlist) > 0) {
            logp <- eval(priorlist$prior.density, as.list(par))
            logpgrad <- unname(sapply(priorlist$prior.grad, function(x, y) ifelse(is.null(x), 0, eval(x, y)), as.list(par)))
        } else {
            logp <- logpgrad <- 0
        }

        if (inherits(v, "try-error")) {
            return(NA)
        } else {
            oldnll <<- v[1] - logp
            grad <- v[-1] * derivpar - logpgrad
            if (length(grad) > 0) names(grad) <- names(derivpar) ## TODO: need a better way of dealing this
            oldgrad <<- grad
            oldpar <<- par

            return(oldnll)
        }
    }

    gradfun <- function(par, data, solver.opts, solver, linklist, priorlist) {
        if (identical(par,oldpar)) {
            return(oldgrad)
        }
        origpar <- apply_link(par, linklist, "linkinv")
        derivpar <- apply_link(par, linklist, "mu.eta")

        v <- try(logLik.sensitivity(origpar, model, data, solver.opts, solver), silent=TRUE)

        if (length(priorlist) > 0) {
            logp <- eval(priorlist$prior.density, as.list(par))
            logpgrad <- unname(sapply(priorlist$prior.grad, function(x, y) ifelse(is.null(x), 0, eval(x, y)), as.list(par)))
        } else {
            logp <- logpgrad <- 0
        }

        if (inherits(v, "try-error")) {
            return(NA)
        } else {
            oldnll <<- v[1] - logp
            grad <- v[-1] * derivpar - logpgrad
            names(grad) <- names(derivpar)
            oldgrad <<- grad
            oldpar <<- par
            return(grad)
        }
    }

    environment(objfun) <- f.env
    environment(gradfun) <- f.env
    parnames <- names(start)
    attr(objfun, "parnames") <- parnames

    if (!keep_sensitivity) gradfun <- NULL

    message("Fitting ode ...")
    m <- mle2(objfun,
              vecpar=TRUE,
              start=start,
              method=method,
              optimizer=optimizer,
              control=control,
              gr=gradfun,
              data=dataarg,
              skip.hessian=skip.hessian,
              ...)

    coef <- apply_link(coef(m), linklist, "linkinv")
    ## calculating hessian in the original scale doesn't make sense when we specify priors
    if (!skip.hessian && length(prior) == 0) {
        if (!length(modelpar)) {
            vcov <- matrix(0, 0, 0)
        } else {
            message("Computing vcov on the original scale ...")

            if (keep_sensitivity && !force.hessian) {
                hessfun <- numDeriv::jacobian
                hmodel <- model
            } else {
                hessfun <- numDeriv::hessian
                hmodel <- Transform(model, keep_sensitivity=FALSE)
            }

            thess <- try(hessfun(logLik.sensitivity, coef,
                                 model=hmodel,
                                 data=data,
                                 solver.opts=solver.opts,
                                 solver=solver,
                                 return.NLL=!keep_sensitivity || force.hessian))
            if(!inherits(thess, "try-error")) {
                if (use.ginv) {
                    vcov <- try(MASS::ginv(thess), silent=TRUE)
                } else {
                    vcov <- try(solve(thess), silent=TRUE)
                }
                if (inherits(vcov, "try-error")) {
                    warning("Couldn't invert Hessian")
                    vcov <- matrix(NA, length(modelpar), length(modelpar))
                }
            } else {
                warning("Couldn't compute Hessian")
                vcov <- matrix(NA,  length(modelpar), length(modelpar))
            }
            rownames(vcov) <- colnames(vcov) <- names(coef)
        }
    } else {
        vcov <- matrix(NA, length(modelpar), length(modelpar))
    }

    new("fitode", model=model, data=data, coef=coef, vcov=vcov,
        min=m@min, mle2=m, link=link,
        fixed=as.list(fixed),
        prior=prior
    )
}

##' Calculate the derivative of an expression with respect to model parameters
##' using sensitivity equations and chain rule
##'
##' @title Calculate the derivative of the mean expression
##' @param model odemodel object
##' @param parms named vector of parameter values
##' @param times time window for which the model should be solved
##' @param solver.opts options for the ode solver (see \code{\link{ode}})
##' @param solver ode solver
ode.sensitivity <- function(model,
                            parms, times,
                            solver.opts=list(method="rk4"),
                            solver=ode) {
    solution <- ode.solve(model, times, parms, solver.opts=solver.opts, solver=solver)

    frame <- c(solution@solution, parms)

    mean <- lapply(model@expr, eval, frame)

    if (model@keep_sensitivity) {
        nstate <- length(model@state)

        sens <- lapply(model@expr.sensitivity, function(expr.sensitivity){
            sens <- matrix(0, nrow=length(times),ncol=length(model@par))

            sens <- Reduce("+", Map("*", lapply(expr.sensitivity$state, eval, frame), solution@sensitivity))

            sens_p <- sapply(expr.sensitivity$par, eval, frame)

            if(is.list(sens_p))
                sens <- sens + do.call("cbind", sens_p)

            sens
        })
    } else {
        sens <-NULL
    }

    list(mean=mean, sensitivity=sens)
}

##' Calculate the derivative of the log-likelihood function with respect to model parameters
##' using sensitivity equations and chain rule
##'
##' @title Calculate the derivative of the log-likelihood function
##' @param parms named vector of parameter values
##' @param model odemodel object
##' @param data data
##' @param solver.opts options for the ode solver (see \code{\link{ode}})
##' @param solver ode solver
##' @param return.NLL (logical) return negative log-likelihood
##' @param return.traj (logical) return estimated trajectory
##' @return a vector of nll and derivative of nll with respect to model parameters
##' (or a list containing (1) the estimated traejctory and (2) a vector of nll and its derivatives)
logLik.sensitivity <- function(parms,
                               model,
                               data,
                               solver.opts=list(method="rk4"),
                               solver=ode,
                               return.NLL=TRUE,
                               return.traj=FALSE) {
    times <- data$times
    ordered.times <- sort(unique(times))

    ss <- ode.sensitivity(model, parms, ordered.times, solver.opts, solver)
    mean <- ss$mean
    sens <- ss$sensitivity

    oo <- match(times, ordered.times)

    nll_list <- sensitivity_list <- vector('list', length(mean))

    for (i in 1:length(nll_list)) {
        ll_fun <- model@loglik[[i]]

        ## skip na observations
        ## this trick allows us to model the difference
        nn <- !is.na(data[,ll_fun@observation])

        frame <- c(list(mean[[i]][oo]), parms, data)

        names(frame)[1] <- ll_fun@mean

        conditional_ll <- eval(ll_fun@expr, frame)[nn]

        nll_list[[i]] <- -sum(conditional_ll)

        if (model@keep_sensitivity) {
            ll_grad <- ll_fun@grad

            loglik.gr <- lapply(ll_grad, eval, frame)

            if (length(model@par) == 1) {
                nll_gr <- -sum((loglik.gr[[1]] * sens[[i]][oo,])[nn])
            } else {
                nll_gr <- -colSums((loglik.gr[[1]] * sens[[i]][oo,])[nn,])
            }

            if(length(loglik.gr) > 1) nll_gr[[ll_fun@par]] <- -sum(loglik.gr[[ll_fun@par]][nn])

            sensitivity_list[[i]] <- nll_gr
        }
    }

    nll <- Reduce("+", nll_list)
    res <- Reduce("+", sensitivity_list)

    if (return.NLL) res <- c(nll, res)

    if (return.traj) res <- list(traj=mean, nll=res)

    res
}
