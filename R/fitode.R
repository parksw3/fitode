##' Set up link functions for ode/loglik parameters
##'
##' @param link named list or vector of strings specifying link functions
##' @param model model.ode object
##' @param loglik loglik.ode object
##' @seealso \code{\link{make.link}}
##' @return list of strings specifying link functions
set_link <- function(link, modelpar) {
    link_default <- as.list(rep("log", length(modelpar)))
    names(link_default) <- modelpar

    if (!missing(link)) link_default[names(link)] <- link

    link_default
}

##' Apply link functions to parameters
##'
##' @param par vector of parameter values
##' @param linklist list containing \code{linkfun}, \code{linkinv}, and \code{mu.eta} for each link
##' @param type string specifying which function should be applied
##' @seealso \code{\link{make.link}}
apply_link <- function(par, linklist, type=c("linkfun", "linkinv", "mu.eta")) {
    type <- match.arg(type)
    ff <- linklist[[type]]
    pp <- unlist(Map(function(x, fun) fun(x), x=par, fun=ff))
    names(pp) <- names(ff)

    pp
}

##' fit ode
##' @rdname fitode
##' @name fitode
##' @param model ode model
##' @param data data frame with time column and observation column
##' @param start named vector of starting parameter values
##' @param tcol time column
##' @param method optimization method
##' @param optimizer optimizer
##' @param link named vector or list of link functions for ode/log-likelihood parameters
##' @param control see optim
##' @param solver.opts options for ode integration. See ode
##' @param skip.hessian skip hessian calculation
##' @param use.ginv use generalized inverse (\code{\link{ginv}}) to compute approximate vcov
##' @param debug print debugging output?
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
                   control=list(maxit=1e5),
                   solver.opts=list(method="rk4"),
                   solver=ode,
                   skip.hessian=FALSE,
                   use.ginv=TRUE,
                   debug=FALSE,
                   ...) {

    modelpar <- model@par

    if ("t" %in% modelpar) {
        stop("`t` is reserved for time variable. Try a different parameterization?")
    }

    if (length(fixed) > 0) {
        fixed <- as.list(fixed)
        if (any(!(names(fixed) %in% modelpar)))
            stop("`fixed`` must be a named vector/list whose names correspond to model parameters")

        tlist <- vector('list', length(fixed))

        for (i in 1:length(fixed)) {
            tlist[[i]] <- as.formula(as.call(c(as.symbol("~"), as.symbol(names(fixed)[i]), unname(fixed[i]))))
        }

        modelpar <- modelpar[!(modelpar %in% names(fixed))]

        model <- Transform(
            model,
            tlist,
            modelpar
        )
    }

    if (!missing(link)) {
        link <- link[!(names(link) %in% names(fixed))]

        if (any(is.na(match(names(link), modelpar)))) stop("Some link functions do not correspond to the model parameters.")
    }

    if (any(is.na(match(modelpar, names(start))))) {

        stop(
            paste0("`start` must specify the following parameters:\n",
                "\node parameters: ", paste(model@par, collapse = ", ")
            )
        )
    }

    ## order parameters ...
    start <- start[modelpar]

    link <- set_link(link, modelpar)

    link_data <- lapply(link, make.link)

    linklist <- lapply(c("linkfun", "linkinv", "mu.eta"),
                       function(x) lapply(link_data, "[[", x))

    names(linklist) <- c("linkfun", "linkinv", "mu.eta")

    newpar <- Map(function(x, y) ifelse(x=="identity", y, paste(x, y, sep=".")), x=link, y=modelpar)
    newpar <- unname(unlist(newpar))

    names(linklist$linkfun) <- names(linklist$mu.eta) <- newpar

    start <- apply_link(start, linklist, "linkfun")

    keep_sensitivity <- model@keep_sensitivity

    names(data)[which(names(data)==tcol)] <- "times"

    dataarg <- list(model=model, data=data, solver.opts=solver.opts, solver=solver, linklist=linklist)

    f.env <- new.env()
    ## set initial values
    assign("oldnll",NULL,f.env)
    assign("oldpar",NULL,f.env)
    assign("oldgrad",NULL,f.env)

    objfun <- function(par, data, solver.opts, solver, linklist) {
        if (identical(par,oldpar)) {
            if (debug) cat("returning old version of value\n")
            return(oldnll)
        }
        if (debug) cat("computing new version (nll)\n")
        origpar <- apply_link(par, linklist, "linkinv")
        derivpar <- apply_link(par, linklist, "mu.eta")

        v <- try(logLik.sensitivity(origpar, model, data, solver.opts, solver), silent=TRUE)
        if (inherits(v, "try-error")) {
            return(NA)
        } else {
            oldnll <<- v[1]
            grad <- v[-1] * derivpar
            if (length(grad) > 0) names(grad) <- names(derivpar) ## TODO: need a better way of dealing this
            oldgrad <<- grad
            oldpar <<- par

            if (debug) {print(oldnll); print(par)}

            return(oldnll)
        }
    }

    gradfun <- function(par, data, solver.opts, solver, linklist) {
        if (identical(par,oldpar)) {
            if (debug) cat("returning old version of grad\n")
            return(oldgrad)
        }
        if (debug) cat("computing new version (grad)\n")
        origpar <- apply_link(par, linklist, "linkinv")
        derivpar <- apply_link(par, linklist, "mu.eta")

        v <- try(logLik.sensitivity(origpar, model, data, solver.opts, solver), silent=TRUE)
        if (inherits(v, "try-error")) {
            return(NA)
        } else {
            oldnll <<- v[1]
            grad <- v[-1] * derivpar
            names(grad) <- names(derivpar)
            oldgrad <<- grad
            oldpar <<- par
            if (debug) {print(oldnll); print(par)}
            return(grad)
        }
    }

    environment(objfun) <- f.env
    environment(gradfun) <- f.env
    parnames <- names(start)
    attr(objfun, "parnames") <- parnames

    if (!keep_sensitivity) gradfun <- NULL ## TODO: I don't like this

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
    if (!skip.hessian && !missing(link)) {
        if (!length(modelpar)) {
            vcov <- matrix(0, 0, 0)
        } else {
            message("Computing vcov on the original scale ...")

            if (keep_sensitivity) {
                hessfun <- numDeriv::jacobian
            } else {
                hessfun <- numDeriv::hessian
            }

            thess <- try(hessfun(logLik.sensitivity, coef,
                                 model=model,
                                 data=data,
                                 solver.opts=solver.opts,
                                 solver=solver,
                                 returnNLL=!keep_sensitivity))
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
        vcov <- vcov(m)
    }

    new("fitode", model=model, data=data, coef=coef, vcov=vcov,
        min=m@min, mle2=m, link=link,
        fixed=fixed
    )
}

##' Calculate sensitivity of the expression with respect to the parameters
##'
##' @param expr expression to be evaluated
##' @param expr.sensitivity partial derivative of expr w.r.t states and parameters
##' @param model model.ode object
##' @param parms named vector of parameter values
##' @param times time window for which the model should be solved
##' @param solver.opts options for the ode solver (see \code{\link{ode}})
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

        ## TODO: implement diffnames
        ## sens <- sens + diff(eval(expr.sensitivity$state[[i]], frame) * solution@sensitivity[[i]])

    } else {
        sens <-NULL
    }

    list(mean=mean, sensitivity=sens)
}

##' Sensitivity of the likelihood function with respect to parameters
##' @param parms named vector of parameter values
##' @param model model.ode object
##' @param data data
##' @param solver.opts options for the ode solver (see \code{\link{ode}})
##' @param returnNLL (logical) return negative log likelihood
##' @return vector of nll and sensitivity of nll with respect to the parameters
logLik.sensitivity <- function(parms,
                               model,
                               data,
                               solver.opts=list(method="rk4"),
                               solver=ode,
                               returnNLL=TRUE) {
    times <- data$times
    ordered.times <- sort(unique(times))

    ss <- ode.sensitivity(model, parms, ordered.times, solver.opts, solver)
    mean <- ss$mean
    sens <- ss$sensitivity

    oo <- match(times, ordered.times)

    nll_list <- sensitivity_list <- vector('list', length(mean))

    for (i in 1:length(nll_list)) {
        ll_fun <- model@loglik[[i]]

        nn <- !is.na(data[,ll_fun@observation])

        frame <- c(list(mean[[i]][oo]), parms, data)

        names(frame)[1] <- ll_fun@mean

        conditional_ll <- eval(ll_fun@expr, frame)[nn]

        nll_list[[i]] <- -sum(conditional_ll)

        if (model@keep_sensitivity) {
            ll_grad <- ll_fun@grad

            loglik.gr <- lapply(ll_grad, eval, frame)

            nll_gr <- -colSums((loglik.gr[[1]] * sens[[i]][oo,])[nn,])

            if(length(loglik.gr) > 1) nll_gr[[ll_fun@par]] <- -sum(loglik.gr[[ll_fun@par]][nn])

            sensitivity_list[[i]] <- nll_gr
        }
    }

    nll <- Reduce("+", nll_list)
    sensitivity <- Reduce("+", sensitivity_list)

    if (!returnNLL) return(sensitivity)

    c(nll, sensitivity)
}

