##' the initializer for model.ode
##'
##' @param .Object object
##' @slot name name of the model
##' @slot model ode model
##' @slot observation observation model
##' @slot initial initial values
##' @slot par parameters
##' @slot keep_sensitivity (logical) maintain the Jacobian as part of the model
##' @examples
##' SI_model <- new("model.ode",
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
##'     par= c("beta", "gamma", "N", "i0", "sigma1", "sigma2")
##' )
##' @docType methods
##' @exportMethod initialize
setMethod(
    "initialize",
    "model.ode",
    function(.Object, name,
             model,
             observation,
             initial,
             par,
             diffnames,
             keep_sensitivity=TRUE) {
        ## TODO: I can't remember why it's asking for a function...
        if (any(sapply(initial, class) != "formula"))
            stop("initial must be a list of formulas or a function")

        .Object@name <- name

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

        if (is.list(observation)) {
            if (any(sapply(observation, class) != "formula"))
                stop("observation must be a list of formulas")

            loglik_list <- lapply(observation, function(ll) {
                ll_model <- select_model(as.character(ll[[3]][[1]]))

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


            .Object@observation <- observation
            .Object@loglik <- lapply(loglik_list, "[[", "ll_model")
            .Object@expr <- lapply(loglik_list, "[[", "expr")
            .Object@expr.sensitivity <- lapply(loglik_list, "[[", "expr.sensitivity")
        } else {
            stop("observation must be a list of formulas or a function")
        }

        if (!missing(diffnames)) .Object@diffnames <- diffnames

        .Object@initial <- initial
        .Object@state <- state
        .Object@par <- par
        .Object@keep_sensitivity <- keep_sensitivity

        .Object
    }
)

##' Evaluate the gradients of a model
##' @param object object to be evaluated
##' @param state state
##' @param par additional parameters
##' @docType methods
##' @exportMethod grad
setMethod(
    "grad",
    "model.ode",
    function(object, state, par) {
        frame <- as.list(c(state, par))
        gr <- sapply(object@grad, eval, frame)
        gr
    }
)

##' Evaluate the jacobian of the gradients
##' @param object object to be evaluated
##' @param state state
##' @param par additional parameters
##' @param type state of par?
##' @docType methods
##' @exportMethod jacobian
setMethod(
    "jacobian",
    "model.ode",
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

##' @exportMethod Transform
setMethod(
    "Transform",
    "model.ode",
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

        new("model.ode",
            object@name,
            newmodel,
            newobservation,
            newinitial,
            par,
            object@diffnames,
            keep_sensitivity)
    }
)

setMethod("show", "model.ode",
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
