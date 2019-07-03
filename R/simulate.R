##' S4 generic for simulating an object
##' @param object an \code{R} object
##' @param ... further arguments passed to methods
setGeneric(
    "simulate",
    def = function(object, ...) {
        standardGeneric("simulate")
    }
)

##' simulate model objects
##'
##' @aliases simulate,fitode-method
##' @param object model.ode object
##' @param times time vector
##' @param parms named vector of parameter values
##' @param y initial values
##' @param solver.opts options for ode solver
##' @param solver ode solver (must take y, times, func, and parms arguments)
##' @param observation (logical) propagate observation error?
##' @import deSolve
##' @docType methods
##' @exportMethod simulate
setMethod("simulate", "model.ode",
    function(object,
             times, parms, y,
             solver.opts=list(method="rk4"),
             solver=ode,
             observation=FALSE,
             nsim=1) {
        simulate_internal(object, times, parms, y, solver.opts, solver, observation, nsim)
    }
)

##' simulate fitode objects
##'
##' @aliases simulate,fitode-method
##' @param object model.ode object
##' @param times time vector
##' @param parms named vector of parameter values
##' @param y initial values
##' @param observation (logical) propagate observation error?
##' @import deSolve
##' @docType methods
##' @exportMethod simulate
setMethod("simulate", "fitode",
    function(object,
             times, parms, y,
             observation=TRUE,
             nsim=1) {
        model <- object@model

        if (missing(parms)) parms <- coef(object)

        if (missing(times)) times <- sort(unique(object@data$times))

        solver.opts <- object@mle2@data$solver.opts

        solver <- object@mle2@data$solver

        simulate_internal(model, times, parms, y, solver.opts, solver, observation, nsim)
    }
)

simulate_internal <- function(object,
                              times, parms, y,
                              solver.opts=list(method="rk4"),
                              solver=ode,
                              observation=TRUE,
                              nsim=1) {
    frame <- as.list(c(parms))

    if (missing(y)) {
        y <- sapply(object@initial, eval, frame)
    } else if (!all(names(y) %in% object@state)) {
        stop("y must have same name as the state variables")
    }

    if (object@keep_sensitivity) {
        ## jacobian(model, state, parms, type="initial")
        ji <- sapply(object@jacobian.initial, function(jj) {
            sapply(jj, eval, frame)
        })
        y <- c(y, ji)
    }

    ss <- new("solution.ode",
              y, times, object, parms,
              solver.opts,
              solver)

    if (!observation) return(ss)

    errorframe <- c(ss@solution, as.list(c(parms)))

    simlist <- vector('list', nsim)

    for (i in 1:nsim) {
        templist <- lapply(object@observation, function(oo) {
            funcall <- oo[[3]]

            erfun <- errorfun(as.character(funcall[[1]]))

            funcall[[1]] <- as.name("erfun")

            eval(funcall, errorframe)
        })

        names(templist) <- sapply(object@observation, function(x) as.character(x[[2]]))
        templist$times <- ss@times
        templist$sim <- i

        simlist[[i]] <- as.data.frame(templist)
    }

    return(do.call("rbind", simlist))
}

errorfun <- function(family=c("dnorm", "dpois", "dnbinom", "dnbinom1", "gamma")) {
    family <- match.arg(family)

    switch(family,
        dnorm=function(mean, sd) rnorm(length(mean), mean=mean, sd=sd),
        dpois=function(lambda) rpois(length(lambda), lambda=lambda),
        dnbinom=function(mu, size) rnbinom(length(mu), mu=mu, size=size),
        dnbinom1=function(mu, phi) rnbinom(length(mu), mu=mu, size=mu/phi),
        dgamma=function(mean, shape) rgamma(length(mean), shape=shape)
    )
}
