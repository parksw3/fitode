##' Class "solution.ode".
##' Result of solving ode modeld with/without sensitivity equations
##'
##' @name solution.ode-class
##' @rdname solution.de-class
##' @slot name name of the model
##' @slot y initial values
##' @slot times time vector
##' @slot model ode model
##' @slot parms parameters of the solution
##' @slot solution solution of the model
##' @slot sensitivity partial derivative of each state variable with respect to the parameters
##' @keywords internal
##' @exportClass solution.ode
setClass(
    "solution.ode",
    slots = c(
        name = "character",
        y = "numeric",
        times = "numeric",
        model = "odemodel",
        parms = "numeric",
        solution = "data.frame",
        sensitivity = "list"
    )
)

##' Constructor for solution.ode class
##'
##' @param .Object object
##' @param y initial values
##' @param times time vector
##' @param model ode model
##' @param parms parameters of the solution
##' @param solver.opts options for ode solver
##' @param solver ode solver (must take y, times, func, and parms as arguments)
##' @return An object of class ``solution.ode'' as described in \code{\link{solution.ode-class}}.
##' @docType methods
##' @keywords internal
##' @exportMethod initialize
setMethod(
    "initialize",
    "solution.ode",
    definition = function(.Object,
                          y, times, model, parms,
                          solver.opts=list(method="rk4"),
                          solver=ode) {
        .Object@name <- model@name
        .Object@y <- y
        .Object@times <- times
        .Object@model <- model
        .Object@parms <- parms

        nstate <- length(model@state)
        npar <- length(model@par)

        diffnames <- model@diffnames

        gfun <- model@gfun

        result <- do.call(solver,
                          c(list(y=y,
                                 times=times,
                                 func=gfun,
                                 parms=parms),
                            solver.opts))

        solution <- as.data.frame(result[,1:(1+nstate)])

        if (length(diffnames) > 0) {
            solution[,diffnames] <- rbind(
                rep(NA, length(diffnames)),
                as.data.frame(diff(as.matrix(solution[,diffnames])))
            )
        }

        .Object@solution <- solution

        if (model@keep_sensitivity) {
            sensitivity <- vector("list", nstate)
            for (i in 1:nstate) {
                sensitivity[[i]] <- result[,(2+nstate):(1+nstate+npar)+(i-1)*npar]
                if(!is.matrix(sensitivity[[i]]))
                    sensitivity[[i]] <- matrix(sensitivity[[i]], ncol=length(model@par))

                if (model@state[i] %in% diffnames)
                    sensitivity[[i]] <- rbind(rep(NA, npar), diff(sensitivity[[i]]))

                colnames(sensitivity[[i]]) <- model@par
            }
            names(sensitivity) <- model@state
            .Object@sensitivity <- sensitivity
        } else {
            .Object@sensitivity <- list()
        }
        .Object
    }
)

##' solve ode models
##' @param model odemodel object
##' @param times time vector
##' @param parms named vector of parameter values
##' @param y initial values
##' @param solver.opts options for ode solver
##' @param solver ode solver (must take y, times, func, and parms as arguments)
##' @return An object of class ``solution.ode'' as described in \code{\link{solution.ode-class}}.
##' @import deSolve
##' @keywords internal
##' @export
ode.solve <- function(model, times, parms, y,
                 solver.opts=list(method="rk4"),
                 solver=ode) {
    frame <- as.list(c(parms))

    if (missing(y)) {
        y <- sapply(model@initial, eval, frame)
    } else if (!all(names(y) %in% model@state)) {
        stop("y must have same name as the state variables")
    }

    if (model@keep_sensitivity) {
        ## jacobian(model, state, parms, type="initial")
        ji <- sapply(model@jacobian.initial, function(jj) {
            sapply(jj, eval, frame)
        })
        y <- c(y, ji)
    }

    new("solution.ode",
        y, times, model, parms,
        solver.opts,
        solver)
}
