##' Set up link functions for ode/loglik parameters
##'
##' @param link named list or vector of strings specifying link functions
##' @param model model.ode object
##' @param loglik loglik.ode object
##' @seealso \code{\link{make.link}}
##' @return list of strings specifying link functions
set_link <- function(link, model, loglik) {
    allpar <- c(model@par, loglik@par)
    link_default <- as.list(rep("identity", length(allpar)))
    names(link_default) <- allpar

    if (!missing(link)) link_default[names(link)] <- link

    loglik.link <- switch(loglik@name,
        gaussian=list(sigma="log"),
        nbinom=list(k="log"),
        nbinom1=list(phi="log")
    )

    if (!is.null(loglik.link)) {
        loglik.link <- loglik.link[!(names(loglik.link) %in% link)]
        link_default[names(loglik.link)] <- loglik.link
    }

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

##' diff function with derivative rules
##' @param x numeric vector
##' @export
.diff <- function(x) diff(x)

##' fit ode
##' @rdname fitode
##' @name fitode
##' @param formula formula specifing observation variable and the mean
##' @param start named vector of starting parameter values
##' @param model ode model
##' @param loglik log liklihood model
##' @param data data frame with time column and observation column
##' @param link named vector or list of link functions for ode/log-likelihood parameters
##' @param control see optim
##' @param ode.opts options for ode integration. See ode
##' @param skip.hessian skip hessian calculation
##' @param use.ginv use generalized inverse (\code{\link{ginv}}) to compute approximate vcov
##' @param debug print debugging output?
##' @param ... mle2 arguments
##' @import bbmle
##' @importFrom numDeriv jacobian hessian
##' @importFrom MASS ginv
##' @seealso \code{\link{mle2}}
##' @export fitode
fitode <- function(formula, start,
                   model, loglik=select_model("gaussian"),
                   data,
                   method="BFGS",
                   optimizer="optim",
                   link,
                   control=list(maxit=1e5),
                   ode.opts=list(method="rk4"),
                   skip.hessian=FALSE,
                   use.ginv=TRUE,
                   debug=FALSE,
                   ...) {

    oldpar <- c(model@par, loglik@par)

    if ("t" %in% oldpar) {
        stop("`t` is reserved for time variable. Try a different parameterization?")
    }

    if (any(is.na(match(names(link), oldpar)))) {
        stop("Some link functions do not correspond to the model parameters.")
    }

    if (any(!is.na(match(loglik@par, model@par)))) {
        stop("Some parameter names in the likeliood model are being used for the model parameters.\nTry a different parameterization?")
    }

    if (any(is.na(match(oldpar, names(start))))) {
        stop(
            paste0("`start` must specify the following parameters:\n",
                "\node parameters: ", paste(model@par, collapse = ", "),
                "\nlikelihood parameters: ", paste(loglik@par)
            )
        )
    }


    ## order parameters ...

    start <- start[oldpar]

    ocol <- as.character(formula[[2]][[2]])
    tcol <- as.character(formula[[2]][[3]])
    data <- data.frame(times = data[[tcol]], observation = data[[ocol]])

    link <- set_link(link, model, loglik)

    link_data <- lapply(link, make.link)

    linklist <- lapply(c("linkfun", "linkinv", "mu.eta"),
                       function(x) lapply(link_data, "[[", x))

    names(linklist) <- c("linkfun", "linkinv", "mu.eta")

    newpar <- Map(function(x, y) ifelse(x=="identity", y, paste(x, y, sep=".")), x=link, y=oldpar)
    newpar <- unname(unlist(newpar))

    names(linklist$linkfun) <- names(linklist$mu.eta) <- newpar

    start <- apply_link(start, linklist, "linkfun")

    keep_sensitivity <- model@keep_sensitivity

    expr <- as.expression(formula[[3]])

    if (keep_sensitivity) {
        expr.sensitivity <- list(
            state=lapply(model@state, function(s) Deriv(expr, s)),
            par=lapply(model@par, function(p) Deriv(expr, p))
        )
    } else {
        expr.sensitivity <- list()
    }

    dataarg <- c(data,list(model=model, loglik=loglik, expr=expr, expr.sensitivity=expr.sensitivity,
                           ode.opts=ode.opts, linklist=linklist, keep_sensitivity=keep_sensitivity))

    ## only accepts one state variable inside .diff
    ## TODO: check that this works...
    if (as.character(expr[[1]][[1]]) == ".diff") {
        dataarg$observation <- dataarg$observation[1:(length(dataarg$observation)-1)]

        if(length(expr[[1]][[2]]) > 1) stop("formula too complicated?")
    }

    f.env <- new.env()
    ## set initial values
    assign("oldnll",NULL,f.env)
    assign("oldpar",NULL,f.env)
    assign("oldgrad",NULL,f.env)
    assign("data", data, f.env)

    objfun <- function(par, expr, expr.sensitivity, model, loglik, observation, times, ode.opts, linklist, keep_sensitivity) {
        if (identical(par,oldpar)) {
            if (debug) cat("returning old version of value\n")
            return(oldnll)
        }
        if (debug) cat("computing new version (nll)\n")
        origpar <- apply_link(par, linklist, "linkinv")
        derivpar <- apply_link(par, linklist, "mu.eta")

        v <- try(logLik.sensitivity(origpar, expr, expr.sensitivity,
                                    model, loglik, observation, times, ode.opts, keep_sensitivity), silent=TRUE)
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

    gradfun <- function(par, expr, expr.sensitivity, model, loglik, observation, times, ode.opts, linklist, keep_sensitivity) {
        if (identical(par,oldpar)) {
            if (debug) cat("returning old version of grad\n")
            return(oldgrad)
        }
        if (debug) cat("computing new version (grad)\n")
        origpar <- apply_link(par, linklist, "linkinv")
        derivpar <- apply_link(par, linklist, "mu.eta")

        v <- try(logLik.sensitivity(origpar, expr, expr.sensitivity,
                                    model, loglik, observation, times, ode.opts, keep_sensitivity), silent=TRUE)
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
        if (!length(oldpar)) {
            vcov <- matrix(0, 0, 0)
        } else {
            message("Computing vcov on the original scale ...")

            if (keep_sensitivity) {
                hessfun <- numDeriv::jacobian
            } else {
                hessfun <- numDeriv::hessian
            }

            thess <- try(hessfun(logLik.sensitivity, coef, expr=expr,
                                            expr.sensitivity=expr.sensitivity,
                                            model=model,loglik=loglik,
                                            observation=dataarg$observation,
                                            times=dataarg$times,
                                            ode.opts=ode.opts,
                                            keep_sensitivity=keep_sensitivity,
                                            returnNLL=!keep_sensitivity))
            if(!inherits(thess, "try-error")) {
                if (use.ginv) {
                    vcov <- try(MASS::ginv(thess), silent=TRUE)
                } else {
                    vcov <- try(solve(thess), silent=TRUE)
                }
                if (inherits(vcov, "try-error")) {
                    warning("Couldn't invert Hessian")
                    vcov <- matrix(NA, length(oldpar), length(oldpar))
                }
            } else {
                warning("Couldn't compute Hessian")
                vcov <- matrix(NA,  length(oldpar), length(oldpar))
            }
            rownames(vcov) <- colnames(vcov) <- names(coef)
        }
    } else {
        vcov <- vcov(m)
    }

    new("fitode", formula=formula, model=model, loglik=loglik,
        data=data, coef=coef, vcov=vcov,
        min=m@min, mle2=m, link=link
    )
}


##' Calculate sensitivity of the expression with respect to the parameters
##'
##' @param expr expression to be evaluated
##' @param expr.sensitivity partial derivative of expr w.r.t states and parameters
##' @param model model.ode object
##' @param parms named vector of parameter values
##' @param times time window for which the model should be solved
##' @param ode.opts options for the ode solver (see \code{\link{ode}})
ode.sensitivity <- function(expr,
                        expr.sensitivity,
                        model,
                        parms, times,
                        ode.opts=list(method="rk4"),
                        keep_sensitivity=TRUE) {
    solution <- ode.solve(model, times, parms, ode.opts=ode.opts)

    frame <- c(solution@solution, parms)

    mean <- eval(expr, frame)

    if (keep_sensitivity) {
        nstate <- length(model@state)

        sens <- matrix(0, nrow=length(mean),ncol=length(model@par))

        if (expr[[1]][[1]]==".diff") {
            for(i in 1:nstate) {
                sens <- sens + diff(eval(expr.sensitivity$state[[i]], frame) * solution@sensitivity[[i]])
            }
        } else {
            for(i in 1:nstate) {
                sens <- sens + eval(expr.sensitivity$state[[i]], frame) * solution@sensitivity[[i]]
            }

            sens_p <- sapply(expr.sensitivity$par, eval, frame)

            if(is.list(sens_p))
                sens <- sens + do.call("cbind", sens_p)
        }
    } else {
        sens=NULL
    }

    list(mean=mean, sensitivity=sens)
}

##' Sensitivity of the likelihood function with respect to parameters
##' @param parms named vector of parameter values
##' @param expr expression to be evaluated
##' @param expr.sensitivity partial derivative of expr w.r.t states and parameters
##' @param model model.ode object
##' @param loglik loglik.ode object
##' @param observation observed values
##' @param times time at which observations were measured
##' @param ode.opts options for the ode solver (see \code{\link{ode}})
##' @param returnNLL (logical) return negative log likelihood
##' @return vector of nll and sensitivity of nll with respect to the parameters
logLik.sensitivity <- function(parms, expr,
                            expr.sensitivity,
                            model, loglik,
                            observation, times=NULL,
                            ode.opts=list(method="rk4"),
                            keep_sensitivity=TRUE,
                            returnNLL=TRUE) {
    if (is.null(times)) times <- seq(length(observation))

    ss <- ode.sensitivity(expr, expr.sensitivity, model, parms, times, ode.opts, keep_sensitivity)
    mean <- ss$mean
    sens <- ss$sensitivity

    loglik.par <- as.list(parms[-c(1:length(model@par))])

    nll <- -sum(Eval(loglik, observation, mean, loglik.par))
    if (keep_sensitivity) {
        loglik.gr <- grad(loglik, observation, mean, loglik.par)
        sensitivity <- c(-colSums(loglik.gr[[1]] * sens))
        if(length(loglik.gr) > 1) sensitivity <- c(sensitivity, -sapply(loglik.gr[-1], sum))
    } else {
        sensitivity <- NULL
    }

    if (!returnNLL) return(sensitivity)

    c(nll, sensitivity)
}
