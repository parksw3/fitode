##' Constructor method of "odemodel" class
##'
##' @param .Object object
##' @param name name of the model
##' @param model ode model
##' @param observation observation model
##' @param initial initial values
##' @param par model parameters
##' @param link link functions for parameters (log links are used as default)
##' @param diffnames optional character vector specifying the names of a variable for which the consecutive difference needs to be calculated
##' @param keep_sensitivity (logical) maintain the Jacobian as a part of the model object?
##' @name odemodel
##' @rdname odemodel-class
##' @examples
##' SI_model <- odemodel(
##'     name = "SI",
##'     model = list(
##'         S ~ - beta*S*I/N,
##'         I ~ beta*S*I/N - gamma*I
##'     ),
##'     observation = list(
##'         susceptible ~ dnorm(mean=S, sd=sigma1),
##'         infected ~ dnorm(mean=I, sd=sigma2)
##'     ),
##'     initial = list(
##'         S ~ N * (1 - i0),
##'         I ~ N * i0
##'     ),
##'     par = c("beta", "gamma", "N", "i0", "sigma1", "sigma2"),
##'     link = c(i0="logit")
##' )
##' @docType methods
##' @exportMethod initialize
setMethod(
    "initialize",
    "odemodel",
    function(.Object, name,
             model,
             observation,
             initial,
             par,
             link,
             diffnames,
             keep_sensitivity=TRUE) {
        ## TODO: I can't remember why it's asking for a function...
        if (any(sapply(initial, class) != "formula"))
            stop("'initial' must be a list of formulas or a function")

        if (any(sapply(observation, class) != "formula"))
            stop("'observation' must be a list of formulas")

        ## warning
        if("dnorm2" %in% sapply(observation, function(ll) as.character(ll[[3]][[1]])) && keep_sensitivity) {
            warning("Sensitivity equations are unavailable for dnorm2 (changing keep_sensitivity=FALSE).")
            keep_sensitivity <- FALSE
        }

        if (missing(name)) name <- "new ODE model"

        .Object@name <- name

        if ("t" %in% par) {
            stop("'t' is reserved for time variable. Try a different parameterization?")
        }

        state <- sapply(initial, function(y) as.character(y[[2]]))
        nstate <- length(state)

        initial <- lapply(initial, function(x) as.expression(x[[3]]))
        names(initial) <- state

        if (keep_sensitivity) {
            deriv <- function(expr, vars) {
                d <- lapply(vars,
                            function(p){
                                Deriv(expr, p)
                            })
                names(d) <- vars
                d
            }

            deriv2 <- function(gradlist, vars) {
                d <- lapply(gradlist,
                            function(x){
                                deriv(x, vars)
                            })
                names(d) <- state
                d
            }
        }

        if (is.list(model)) {
            if (any(sapply(model, class) != "formula"))
                stop("model must be a list of formulas or a function")

            grad <- lapply(model, function(x) as.expression(x[[3]]))
            names(grad) <- state

            if (keep_sensitivity) {
                .Object@jacobian.initial <- deriv2(initial, par)
                .Object@jacobian.state <- jacobian.state <- deriv2(grad, state)
                .Object@jacobian.par <- jacobian.par <- deriv2(grad, par)
            }

            if(keep_sensitivity) {
                gfun <- function(times, y, parms) {
                    state <- y[1:nstate]
                    frame <- as.list(c(t=times, state, parms))
                    ## equivalent to `grad(model, state, parms)` but faster
                    gr <- sapply(grad, eval, frame)
                    ## jacobian(model, state, parms, type="state")
                    js <- sapply(jacobian.state, function(jj) {
                        sapply(jj, eval, frame)
                    })
                    ## jacobian(model, state, parms, type="par")
                    jp <- sapply(jacobian.par, function(jj) {
                        sapply(jj, eval, frame)
                    })
                        list(c(gr, matrix(y[-c(1:nstate)], ncol=nstate) %*% js + jp))
                }
            } else {
                gfun <- function(times, y, parms) {
                    frame <- as.list(c(t=times, y, parms))
                    gr <- sapply(grad, eval, frame)
                        list(c(gr))
                }
            }

            .Object@grad <- grad
            .Object@gfun <- gfun

        } else if (is.function(model)) {
            keep_sensitivity <- FALSE
            .Object@gfun <- model
        } else {
            stop("model must be a list of formulas or a function")
        }

        loglik_list <- lapply(observation, function(ll) {
            if (as.character(ll[[3]][[1]])=="ols" && length(observation) > 1) {
                stop("'ols' is only available for univariate time series")
            }

            ll_model <- select_model(as.character(ll[[3]][[1]]))

            ll_check <- c(match(c(ll_model@mean, ll_model@par), names(as.list(ll[[3]])[-1])),
                          match(names(as.list(ll[[3]])[-1]), c(ll_model@mean, ll_model@par)))

            if (any(is.na(ll_check))) {
                ll_check_msg <- paste0(
                    "'", ll_model@name, "' requires following arguments:\n",
                    "'", ll_model@mean, "' for specifying the mean trajectory"
                )

                if (length(ll_model@par) > 0) {
                    ll_check_msg <- paste0(
                        ll_check_msg,
                        " and '", ll_model@par, "' for specifying the amount of dispersion"
                    )
                }

                stop(ll_check_msg)
            }

            trans_obs <- as.formula(as.call(c(as.symbol("~"), as.symbol("X"), ll[[2]])))

            trans_list <- list(trans_obs)

            likpar <- ll_model@par

            if (length(likpar) > 0) {
                call <- as.list(ll[[3]])[[likpar]]

                trans_list <- append(trans_list, as.formula(as.call(c(as.symbol("~"), as.symbol(likpar), call))))
            } else {
                call <- likpar
            }

            ll_model <- Transform(ll_model,
                                  observation=as.character(ll[[2]]),
                                  transforms=trans_list,
                                  par=call,
                                  keep_grad=keep_sensitivity
            )

            expr <- ll[[3]][[ll_model@mean]]

            if (keep_sensitivity) {
                expr.sensitivity <- list(
                    state=lapply(state, function(s) Deriv(expr, s)),
                    par=lapply(par, function(p) Deriv(expr, p))
                )

                names(expr.sensitivity$state) <- state
                names(expr.sensitivity$par) <- par
            } else {
                expr.sensitivity <- list()
            }

            list(ll_model=ll_model,
                 expr=expr,
                 expr.sensitivity=expr.sensitivity)
        })

        ## set up link functions
        if (!missing(link)) {
            link <- link[names(link) %in% par]

            if (any(is.na(match(names(link), par)))) stop("Some link functions do not correspond to the model parameters.")
        }

        link <- unlist(set_link(link, par))
        .Object@link <- link

        .Object@observation <- observation
        .Object@loglik <- lapply(loglik_list, "[[", "ll_model")
        .Object@expr <- lapply(loglik_list, "[[", "expr")
        .Object@expr.sensitivity <- lapply(loglik_list, "[[", "expr.sensitivity")

        if (!missing(diffnames)) .Object@diffnames <- diffnames

        .Object@initial <- initial
        .Object@state <- state
        .Object@par <- par
        .Object@keep_sensitivity <- keep_sensitivity

        .Object
    }
)

##' Wrapper function odemodel
##' @name odemodel
##' @rdname odemodel-class
##' @keywords internal
##' @export
odemodel <- function(...) new("odemodel", ...)

##' Evaluate the gradients of a model
##' @param object odemodel object
##' @param state state
##' @param par parameter values
##' @docType methods
##' @keywords internal
##' @exportMethod grad
setMethod(
    "grad",
    "odemodel",
    function(object, state, par) {
        frame <- as.list(c(state, par))
        gr <- sapply(object@grad, eval, frame)
        gr
    }
)

##' Evaluate the jacobian of the gradients
##' @param object odemodel object
##' @param state state
##' @param par parameter values
##' @param type state of par?
##' @docType methods
##' @keywords internal
##' @exportMethod jacobian
setMethod(
    "jacobian",
    "odemodel",
    definition <- function(object, state, par, type=c("initial", "state", "par")) {
        type <- match.arg(type)
        frame <- as.list(c(state, par))
        jc <- switch(type,
            initial=object@jacobian.initial,
            state=object@jacobian.state,
            par=object@jacobian.par
        )

        l <- sapply(jc, function(jj) {
            sapply(jj, eval, frame)
        })
        l
    }
)

##' Transform the model
##' @param object odemodel object
##' @param transforms list of formulas specifying transformations
##' @param par model parameters
##' @param keep_sensitivity (logical) maintain the Jacobian as part of the model
##' @keywords internal
##' @exportMethod Transform
setMethod(
    "Transform",
    "odemodel",
    function(object, transforms=NULL, par, keep_sensitivity) {

        if (missing(keep_sensitivity)) keep_sensitivity <- object@keep_sensitivity

        if (is.null(transforms)) transforms <- list()

        allvars <- c(object@par)
        transforms <- trans(transforms, allvars)

        nstate <- length(object@state)

        newmodel <- newinitial <- vector('list', nstate)

        if (length(object@grad) > 0) {
            for(i in 1:nstate) {
                fixed <- c(as.symbol("~"), as.symbol(object@state[i]))
                ff <- lapply(list(object@grad[[i]], object@initial[[i]]), function(x){
                    f <- c(fixed, subst(x[[1]], transforms))
                    f <- as.formula(as.call(f))
                })
                newmodel[[i]] <- ff[[1]]
                newinitial[[i]] <- ff[[2]]
            }
        }

        newobservation <- lapply(object@observation, function(x) {
            x[[3]] <- as.call(lapply(as.list(x[[3]]), subst, transforms))
            x
        })

            ## TODO: what happens when length(object@grad) > 0??

        if (missing(par)) par <- object@par

        new("odemodel",
            object@name,
            newmodel,
            newobservation,
            newinitial,
            par,
            object@link,
            object@diffnames,
            keep_sensitivity)
    }
)

##' Show the model
##' @param object odemodel object
##' @keywords internal
##' @exportMethod show
setMethod("show", "odemodel",
    function(object){
        cat("Name:", object@name, "\n")

        cat("\nObservations:\n")
        for(i in 1:length(object@observation)) {
            cat(deparse(object@observation[[i]]), "\n")
        }

        cat("\nInitial values:\n")
        g <- paste0(object@state, "(0) = ", sapply(object@initial, as.character))
        for(i in 1:length(g))
            cat(g[i], "\n")

        cat("\nParameters:", object@par, "\n")
    }
)
