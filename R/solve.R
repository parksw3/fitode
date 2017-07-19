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
                sensitivity[[i]] <- result[,(2+nstate):(1+nstate+npar)+(i-1)*npar]
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

##' @import deSolve
##' @export
solve <- function(model, times, parms,
                 y,
                 keep_sensitivity=TRUE,
                 test=FALSE,
                 ...) {
    if (missing(y)) {
        frame <- as.list(c(parms))
        y <- sapply(model@initial, eval, frame)
    } else if (names(y) != model@state) {
        stop("y must have same name as the state variables")
    }

    if (keep_sensitivity) {
        nstate <- length(model@state)

        jacobian(model, y, parms, "initial")

        ## jacobian(model, state, parms, type="initial")
        ji <- sapply(model@jacobian.initial, function(jj) {
            sapply(jj, eval, frame)
        })
        yini <- c(y, ji)

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
        yini <- y
        gfun <- function(y, times, parms) {
            gr <- grad(model, y, parms)
        }
    }

    new("solution.ode",
        model,
        ode(yini, times, gfun, parms, method="rk4", hini=0.1),
        keep_sensitivity)
}
