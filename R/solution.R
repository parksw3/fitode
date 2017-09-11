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
##' @slot keep_sensitivity keep sensitivity equations
##' @exportClass solution.ode
setClass(
    "solution.ode",
    slots = c(
        name = "character",
        y = "numeric",
        times = "numeric",
        model = "model.ode",
        parms = "numeric",
        solution = "data.frame",
        sensitivity = "list",
        keep_sensitivity ="logical"
    )
)

##' the initializer for model.ode
##'
##' @param .Object object
##' @param y initial values
##' @param times time vector
##' @param model ode model
##' @param parms parameters of the solution
##' @param ode.opts options for ode integration
##' @param keep_sensitivity keep sensitivity equations
##' @docType methods
##' @exportMethod initialize
setMethod(
    "initialize",
    "solution.ode",
    definition = function(.Object,
                          y, times, model, parms,
                          ode.opts=list(method="rk4", hini=0.1),
                          keep_sensitivity=TRUE) {
        .Object@name <- model@name
        .Object@y <- y
        .Object@times <- times
        .Object@model <- model
        .Object@parms <- parms

        nstate <- length(model@state)
        npar <- length(model@par)

        if (keep_sensitivity) {
            gfun <- function(times, y, parms) {
                state <- y[1:nstate]
                frame <- as.list(c(state, parms))
                ## equivalent to `grad(model, state, parms)` but faster
                gr <- sapply(model@grad, eval, frame)
                ## jacobian(model, state, parms, type="state")
                js <- sapply(model@jacobian.state, function(jj) {
                    sapply(jj, eval, frame)
                })
                ## jacobian(model, state, parms, type="par")
                jp <- sapply(model@jacobian.par, function(jj) {
                    sapply(jj, eval, frame)
                })

                list(c(gr, matrix(y[-c(1:nstate)], ncol=nstate) %*% js + jp))
            }
        } else {
            gfun <- function(times, y, parms) {
                frame <- as.list(c(y, parms))
                gr <- sapply(model@grad, eval, frame)

                list(c(gr))
            }
        }

        result <- do.call("ode",
                          c(list(y=y,
                                 times=times,
                                 func=gfun,
                                 parms=parms),
                            ode.opts))

            ode(y, times, gfun, parms, method="rk4", hini=0.1)

        .Object@solution <- as.data.frame(result[,1:(1+nstate)])

        if (keep_sensitivity) {
            sensitivity <- vector("list", nstate)
            for (i in 1:nstate) {
                sensitivity[[i]] <- result[,(2+nstate):(1+nstate+npar)+(i-1)*npar]
                if(!is.matrix(sensitivity[[i]]))
                    sensitivity[[i]] <- matrix(sensitivity[[i]], ncol=length(model@par))

                colnames(sensitivity[[i]]) <- model@par
            }
            names(sensitivity) <- model@state
            .Object@sensitivity <- sensitivity
        } else {
            .Object@sensitivity <- list()
        }
        .Object@keep_sensitivity <- keep_sensitivity
        .Object
    }
)

##' solve ode models
##' @param model model.ode object
##' @param times time vector
##' @param parms named vector of parameter values
##' @param y initial values
##' @param ode.opts options for ode integration
##' @param keep_sensitivity keep sensitivity equations
##' @import deSolve
##' @export
ode.solve <- function(model, times, parms, y,
                 ode.opts=list(method="rk4", hini=0.1),
                 keep_sensitivity=TRUE) {
    if (missing(y)) {
        frame <- as.list(c(parms))
        y <- sapply(model@initial, eval, frame)
    } else if (names(y) != model@state) {
        stop("y must have same name as the state variables")
    }

    if (keep_sensitivity) {
        nstate <- length(model@state)

        ## jacobian(model, state, parms, type="initial")
        ji <- sapply(model@jacobian.initial, function(jj) {
            sapply(jj, eval, frame)
        })
        y <- c(y, ji)
    }

    new("solution.ode",
        y, times, model, parms,
        ode.opts,
        keep_sensitivity)
}
