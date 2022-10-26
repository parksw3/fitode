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
##' @param call original function call
##' @return An object of class ``odemodel'' as described in \code{\link{odemodel-class}}.
## removed @name, @rdname as otherwise the appropriate \alias{} doesn't show up and we get a warning about an undocumented method ...
## https://stackoverflow.com/questions/7356120/how-to-properly-document-s4-methods-using-roxygen2
## says we should need only @rdname and @aliases ...
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
##' @importFrom methods initialize
##' @export initialize
##' @export
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
             keep_sensitivity=TRUE,
             call) {
        ## TODO: I can't remember why it's asking for a function...
        if (any(sapply(initial, class) != "formula"))
            stop("'initial' must be a list of formulas or a function")

        if (any(sapply(observation, class) != "formula"))
            stop("'observation' must be a list of formulas")

        ## warning
        if("dnorm2" %in% vapply(observation, get_head, character(1)) && keep_sensitivity) {
            warning("Sensitivity equations are unavailable for dnorm2 (setting keep_sensitivity=FALSE).")
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

        models <- vapply(observation, get_head, character(1))
        if (any(models=="ols") && !all(models=="ols")) {
            stop("ols() must apply to all observations within a model")
        }
        lfun <- function(ll) {

            ll_model <- select_model(get_head(ll))
            ## FIXME: allow matching by position??
            rhs_names <- names(as.list(get_rhs(ll))[-1])
            ll_check <- c(match(c(ll_model@mean, ll_model@par), rhs_names),
                          match(rhs_names, c(ll_model@mean, ll_model@par)))

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
                Lcall <- as.list(ll[[3]])[[likpar]]
                trans_list <- append(trans_list, as.formula(as.call(c(as.symbol("~"), as.symbol(likpar), Lcall))))
            } else {
                Lcall <- likpar
            }
            ll_model <- Transform(ll_model,
                                  observation=deparse(ll[[2]]),
                                  transforms=trans_list,
                                  par=Lcall,
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
        } ## lfun
        loglik_list <- lapply(observation, lfun)

        ## set up link functions
        ## FIXME:: how does this compare to check_link?
        if (!missing(link)) {

            nomatch <- which(is.na(match(names(link), par)))
            if (length(nomatch)>0) stop("Some link functions do not correspond to the model parameters: ",
                                        paste(names(link)[nomatch],collapse=", "))

            link <- link[names(link) %in% par]

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
        .Object@call <- call

        .Object
    }
)

##' Create a new odemodel
##' @name odemodel
##' @rdname odemodel-class
##' @keywords internal
##' @export
odemodel <- function(...) {
    call <- match.call()
    new("odemodel", call=call, ...)
}

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
##' @param observation observation model
##' @param initial initial values
##' @param par model parameters
##' @param keep_sensitivity (logical) maintain the Jacobian as part of the model
##' @return An object of class ``odemodel'' as described in \code{\link{odemodel-class}}.
##' @keywords internal
##' @exportMethod Transform
setMethod(
    "Transform",
    "odemodel",
    function(object, transforms, observation, initial, par, link, keep_sensitivity) {
        if (missing(keep_sensitivity)) keep_sensitivity <- object@keep_sensitivity

        if (missing(transforms)) transforms <- list()

        if (missing(observation)) observation <- object@observation

        if (missing(initial)) initial <- object@initial

        if (missing(link)) link <- object@link

        allvars <- c(object@par)
        transforms <- trans(transforms, allvars)

        nstate <- length(object@state)

        newmodel <- newinitial <- vector('list', nstate)

        if (length(object@grad) > 0) {
            for(i in 1:nstate) {
                fixed <- c(as.symbol("~"), as.symbol(object@state[i]))
                ff <- lapply(list(object@grad[[i]], initial[[i]]), function(x){
                    f <- c(fixed, subst(x[[1]], transforms))
                    f <- as.formula(as.call(f))
                })
                newmodel[[i]] <- ff[[1]]
                newinitial[[i]] <- ff[[2]]
            }
        }

        newobservation <- lapply(observation, function(x) {
            x[[3]] <- as.call(lapply(as.list(x[[3]]), subst, transforms))
            x
        })

        if (missing(par)) par <- object@par

        odemodel(
            object@name,
            newmodel,
            newobservation,
            newinitial,
            par,
            link,
            object@diffnames,
            keep_sensitivity)
    }
)

##' Show the model
##' @param object odemodel object
##' @return No return value, called for side effects
##' @keywords internal
##' @importFrom methods show
##' @export show
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
