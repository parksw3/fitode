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

##' fit ode
##' @rdname fitode
##' @name fitode
##' @param formula formula specifing observation variable and the mean
##' @param start named vector of starting parameter values
##' @param model ode model
##' @param loglik log liklihood model
##' @param data data frame with time column and observation column
##' @param tcol time column
##' @param link named vector or list of link functions for ode/log-likelihood parameters
##' @param control see optim
##' @param ode.opts options for ode integration. See ode
##' @param debug print debugging output?
##' @import bbmle
##' @importFrom numDeriv jacobian
##' @export fitode
setMethod(
    "initialize",
    "fitode",
    function(.Object,
             formula, start,
             model, loglik=select_model("gaussian"),
             data,
             method="BFGS",
             optimizer="optim",
             tcol = "times",
             link,
             control=list(maxit=1e5),
             ode.opts=list(method="lsoda"),
             skip.hessian=FALSE,
             debug=FALSE) {
        oldpar <- c(model@par, loglik@par)

        if (any(is.na(match(names(start), oldpar)))) {
            stop(
                paste0("`start` must specify the following parameters:\n",
                    "\node parameters: ", paste(model@par, collapse = ", "),
                    "\nlikelihood parameters: ", paste(loglik@par)
                )
            )
        }

        ## order parameters ...

        start <- start[oldpar]

        .Object@model <- model
        .Object@formula <- formula
        .Object@loglik <- loglik

        ocol <- as.character(formula[[2]])
        data <- data.frame(times = data[[tcol]], observation = data[[ocol]])

        .Object@data <- data

        link <- set_link(link, model, loglik)

        link_data <- lapply(link, make.link)

        linklist <- lapply(c("linkfun", "linkinv", "mu.eta"),
                           function(x) lapply(link_data, "[[", x))

        names(linklist) <- c("linkfun", "linkinv", "mu.eta")

        .Object@link <- link

        newpar <- Map(function(x, y) ifelse(x=="identity", y, paste(x, y, sep=".")), x=link, y=oldpar)
        newpar <- unname(unlist(newpar))

        names(linklist$linkfun) <- names(linklist$mu.eta) <- newpar
        start <- apply_link(start, linklist, "linkfun")

        dataarg <- c(data,list(model=model, loglik=loglik, formula=formula, ode.opts=ode.opts, linklist=linklist))

        f.env <- new.env()
        ## set initial values
        assign("oldnll",NULL,f.env)
        assign("oldpar",NULL,f.env)
        assign("oldgrad",NULL,f.env)
        assign("data", data, f.env)

        objfun <- function(par, formula, model, loglik, observation, times, ode.opts, linklist) {
            if (identical(par,oldpar)) {
                if (debug) cat("returning old version of value\n")
                return(oldnll)
            }
            if (debug) cat("computing new version (nll)\n")
            origpar <- apply_link(par, linklist, "linkinv")
            derivpar <- apply_link(par, linklist, "mu.eta")

            v <- try(logLik.sensitivity(origpar, formula, model, loglik, observation, times, ode.opts))
            if (inherits(v, "try-error")) {
                return(NA)
            } else {
                oldnll <<- v[1]
                grad <- v[-1] * derivpar
                names(grad) <- names(derivpar)
                oldgrad <<- grad
                oldpar <<- par

                return(oldnll)
            }
        }
        gradfun <- function(par, formula, model, loglik, observation, times, ode.opts, linklist) {
            if (identical(par,oldpar)) {
                if (debug) cat("returning old version of grad\n")
                return(oldgrad)
            }
            if (debug) cat("computing new version (grad)\n")
            origpar <- apply_link(par, linklist, "linkinv")
            derivpar <- apply_link(par, linklist, "mu.eta")

            v <- try(logLik.sensitivity(origpar, formula, model, loglik, observation, times, ode.opts))
            if (inherits(v, "try-error")) {
                return(NA)
            } else {
                oldnll <<- v[1]
                grad <- v[-1] * derivpar
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

        message("Fitting ode ...")

        m <- mle2(objfun,
                  vecpar=TRUE,
                  start=start,
                  method=method,
                  optimizer=optimizer,
                  control=control,
                  gr=gradfun,
                  data=dataarg,
                  skip.hessian=skip.hessian)

        .Object@mle2 <- m

        coef <- apply_link(coef(m), linklist, "linkinv")

        if (!skip.hessian && !missing(link)) {
            message("Computing vcov on the original scale ...")
            thess <- numDeriv::jacobian(logLik.sensitivity, coef, formula=formula,
                model=model,loglik=loglik,
                observation=data[,2],
                times=data[,1],
                ode.opts=ode.opts,
                returnNLL=FALSE)
            vcov <- solve(thess)
            colnames(vcov) <- rownames(vcov) <- names(coef)
        } else {
            vcov <- vcov(m)
        }

        .Object@coef <- coef
        .Object@vcov <- vcov
        .Object@min <- m@min

        .Object
    }
)

##' Calculate sensitivity of the expression with respect to the parameters
##'
##' @param expr expression of model states and parameters
##' @param model model.ode object
##' @param parms named vector of parameter values
##' @param times time window for which the model should be solved
##' @param ode.opts options for the ode solver (see \code{\link{ode}})
##' @examples
##' SI_model <- new("model.ode",
##'     name = "SI",
##'     model = list(
##'         S ~ - beta*S*I/N,
##'         I ~ beta*S*I/N - gamma*I
##'     ),
##'     initial = list(
##'         S ~ N * (1 - i0),
##'         I ~ N * i0
##'     ),
##'     par= c("beta", "gamma", "N", "i0")
##' )
##'
##' ode.sensitivity(expression(gamma*I), SI_model, parms=c(beta=2, gamma=1, N=1e5, i0=1e-4), times=1:10)
##' @export
ode.sensitivity <- function(expr, model,
                        parms, times,
                        ode.opts=list(method="lsoda")) {
    capture.output(solution <- ode.solve(model, times, parms, ode.opts=ode.opts))

    frame <- c(solution@solution, parms)

    mean <- eval(expr, frame)

    nstate <- length(model@state)

    dmds <- lapply(model@state, function(s) Deriv(expr, s))
    dmdp <- lapply(model@par, function(p) Deriv(expr, p))

    sens <- matrix(0, nrow=length(times),ncol=length(model@par))
    for(i in 1:nstate) {
        sens <- sens + eval(dmds[[i]], frame) * solution@sensitivity[[i]]
    }

    sens_p <- sapply(dmdp, eval, frame)

    if(is.list(sens_p))
        sens <- sens + do.call("cbind", sens_p)

    list(mean=mean, sensitivity=sens)
}

##' Sensitivity of the likelihood function with respect to parameters
##' @param parms named vector of parameter values
##' @param formula formula specifing observation variable and the mean
##' @param model model.ode object
##' @param loglik loglik.ode object
##' @param observation observed values
##' @param times time at which observations were measured
##' @param ode.opts options for the ode solver (see \code{\link{ode}})
##' @param returnNLL (logical) return negative log likelihood
##' @return vector of nll and sensitivity of nll with respect to the parameters
logLik.sensitivity <- function(parms, formula,
                            model, loglik,
                            observation, times=NULL,
                            ode.opts=list(method="lsoda"),
                            returnNLL=TRUE) {
    if (is.null(times)) times <- seq(length(observation))
    expr <- as.expression(formula[[3]])
    ss <- ode.sensitivity(expr, model, parms, times, ode.opts)
    mean <- ss$mean
    sens <- ss$sensitivity

    loglik.par <- as.list(parms[-c(1:length(model@par))])

    nll <- -sum(Eval(loglik, observation, mean, loglik.par))
    loglik.gr <- grad(loglik, observation, mean, loglik.par)
    sensitivity <- c(-colSums(loglik.gr[[1]] * sens))
    if(length(loglik.gr) > 1) sensitivity <- c(sensitivity, -sapply(loglik.gr[-1], sum))

    if (!returnNLL) return(sensitivity)

    c(nll, sensitivity)
}
