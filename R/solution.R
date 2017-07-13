setClass(
    "solution.ode",
    slots = c(
        name = "character",
        model = "model.ode",
        solution = "data.frame",
        sensitivity = "list",
        keep_sensitivity ="logical"
    )
)

##' the initializer for model.ode
##'
##' @param .Object object
##' @slot model ode model
##' @slot result deSolve result
##' @docType methods
##' @exportMethod initialize
setMethod(
    "initialize",
    "solution.ode",
    definition = function(.Object,
                          model,
                          result,
                          keep_sensitivity=TRUE) {
        .Object@name <- model@name
        .Object@model <- model
        nstate <- length(model@state)
        npar <- length(model@par)
        .Object@solution <- as.data.frame(result[,1:(1+nstate)])

        if (keep_sensitivity) {
            sensitivity <- vector("list", nstate)
            for (i in 1:nstate) {
                sensitivity[[i]] <- as.data.frame(result[,(2+nstate):(1+nstate+npar)+(i-1)*npar])
                names(sensitivity[[i]]) <- model@par
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

##' @import deSolve
##' @export
solve <- function(model, times, parms,
                 y,
                 method="lsoda",
                 keep_sensitivity=TRUE,
                 ...) {
    if (missing(y)) {
        frame <- as.list(c(parms))
        y <- sapply(model@initial, eval, frame)
    } else if (names(y) != model@state) {
        stop("y must have same name as the state variables")
    }

    if (keep_sensitivity) {
        nstate <- length(model@state)
        nsens <- nstate * length(model@par)

        jacobian(model, y, parms, "initial")

        ## TODO: allow for expressions in y so that it can be affected by the parameters
        yini <- c(y, jacobian(model, y, parms, type="initial"))

        gfun <- function(times, y, parms) {
            state <- y[1:nstate]
            gr <- grad(model, state, parms)
            js <- jacobian(model, state, parms, type="state")
            jp <- jacobian(model, state, parms, type="par")

            list(c(gr, matrix(y[-c(1:nstate)], ncol=nstate) %*% js + jp))
        }

    } else {
        yini <- y
        gfun <- function(y, times, parms) {
            gr <- grad(model, y, parms)
        }
    }

    new("solution.ode",
        model,
        ode(yini, times, gfun, parms, method),
        keep_sensitivity)
}

