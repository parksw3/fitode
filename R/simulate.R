##' simulate model objects
##'
##' @aliases simulate,odemodel-method
##' @param object odemodel object
##' @param nsim number of simulations
##' @param seed random-number seed
##' @param times time vector
##' @param parms named vector of parameter values
##' @param y initial values
##' @param solver.opts options for ode solver
##' @param solver ode solver (must take y, times, func, and parms arguments)
##' @param observation (logical) propagate observation error?
##' @return The numerical simulation of the object
##' @docType methods
##' @importFrom bbmle simulate
##' @exportMethod simulate
setMethod("simulate", "odemodel",
    function(object, nsim=1,
             seed=NULL,
             times, parms, y,
             solver.opts=list(method="rk4"),
             solver=ode,
             observation=TRUE) {
        simulate_internal(object, times, parms, y, solver.opts, solver, observation, nsim, seed)
    }
)

##' simulate fitode objects
##'
## FIXME: use @inheritParams ?
##' @aliases simulate,fitode-method
##' @param object fitode object
##' @param nsim number of simulations
##' @param seed random-number seed
##' @param times time vector
##' @param parms named vector of parameter values
##' @param y initial values
##' @param observation (logical) propagate observation error?
##' @return The numerical simulation of the object
##' @docType methods
setMethod("simulate", "fitode",
    function(object, nsim=1,
             seed=NULL,
             times, parms, y,
             observation=TRUE) {
        model <- object@model

        if (missing(parms)) parms <- coef(object)

        if (missing(times)) times <- sort(unique(object@data$times))

        solver.opts <- object@mle2@data$solver.opts

        solver <- object@mle2@data$solver

        simulate_internal(model, times, parms, y, solver.opts, solver, observation, nsim, seed)
    }
)

##' Simulates deterministic trajectories and propagates observation error
##'
##' @title Internal function for simulation models
##' @param model odemodel object
##' @param times time vector
##' @param parms named vector of parameter values
##' @param y initial values
##' @param solver.opts options for ode solver
##' @param solver ode solver (must take y, times, func, and parms arguments)
##' @param observation (logical) propagate observation error?
##' @param nsim number of simulations
##' @param seed seed
simulate_internal <- function(model,
                              times, parms, y,
                              solver.opts=list(method="rk4"),
                              solver=ode,
                              observation=TRUE,
                              nsim=1,
                              seed=NULL) {
    if (!is.null(seed)) set.seed(seed)

    ss <- ode.solve(model, times, parms, y, solver.opts, solver)

    out <- as.data.frame(ss@solution)

    if (!observation) return(out)

    errorframe <- c(out, as.list(c(parms)))

    simlist <- vector('list', nsim)

    for (i in 1:nsim) {
        templist <- lapply(model@observation, function(oo) {
            funcall <- oo[[3]]

            erfun <- errorfun(as.character(funcall[[1]]))

            funcall[[1]] <- as.name("erfun")

            eval(funcall, errorframe)
        })

        names(templist) <- sapply(model@observation, function(x) as.character(x[[2]]))
        templist$time <- out$time
        templist$sim <- i

        simlist[[i]] <- as.data.frame(templist)
    }

    out2 <- do.call("rbind", simlist)

    out3 <- merge(out, out2, by="time")
    out3 <- out3[order(out3$sim),]

    names(out3)[1] <- "times"

    return(out3)
}

#' @importFrom stats rnorm rpois rnbinom rgamma rlnorm
errorfun <- function(family=c("ols", "dnorm", "dpois", "dnbinom", "dnbinom1", "dgamma", "dlnorm")) {
    family <- match.arg(family)

    switch(family,
        ols=function(mean) rnorm(length(mean), mean=mean, sd=1), ## FIXME: should this be zero-sd? or empirical sd?
        dnorm=function(mean, sd) rnorm(length(mean), mean=mean, sd=sd),
        dpois=function(lambda) rpois(length(lambda), lambda=lambda),
        dnbinom=function(mu, size) rnbinom(length(mu), mu=mu, size=size),
        dnbinom1=function(mu, phi) rnbinom(length(mu), mu=mu, size=mu/phi),
        dgamma=function(mean, shape) rgamma(length(mean), shape=shape),
        dlnorm=function(meanlog, sdlog) rlnorm(length(meanlog), meanlog=meanlog,
                                               sdlog=sdlog)
    )
}
