##' Class representing ode models
##'
##' @slot name name of the model
##' @slot grad list of expressions representing the gradients
##' @slot initial list of expressions representing the initial values
##' @slot state state variables
##' @slot par parameters
##' @slot jacobian.initial Jacobian of initial values with respect to its parameters
##' @slot jacobian.state Jacobian with respect to its states
##' @slot jacobian.par Jacobian with repsect to its parameters
setClass(
    "model.ode",
    slots = c(
        name = "character",
        grad = "list",
        initial= "list",
        state = "character",
        par = "character",
        jacobian.initial = "list",
        jacobian.state = "list",
        jacobian.par = "list"
    )
)

##' the initializer for model.ode
##'
##' @param .Object object
##' @slot name name of the model
##' @slot model ode model
##' @slot initial initial values
##' @slot par parameters
##' @slot keep_jacobian (logical) maintain the Jacobian as part of the model
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
##' @docType methods
##' @exportMethod initialize
setMethod(
    "initialize",
    "model.ode",
    definition = function(.Object, name,
                          model, initial,
                          par,
                          keep_jacobian=TRUE) {
        .Object@name <- name

        if (any(sapply(model, class) != "formula"))
            stop("model must be a list of formulas")

        if (any(sapply(initial, class) != "formula"))
            stop("initial must be a list of formulas")

        mi <- list(model, initial)

        state <- lapply(mi, function(x){
            sapply(x, function(y) as.character(y[[2]]))
        })

        if (!do.call(all.equal, state)) {
            stop("initial values do not have same state variable names as the model provided")
        } else {
            state <- state[[1]]
        }


        grad <- lapply(model, function(x) as.expression(x[[3]]))
        initial <- lapply(initial, function(x) as.expression(x[[3]]))
        names(grad) <- state
        names(initial) <- state

        .Object@grad <- grad
        .Object@initial <- initial
        .Object@state <- state
        .Object@par <- par

        deriv <- function(expr, vars) {
            d <- lapply(vars,
                        function(p){
                            D(expr, p)
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

        if (keep_jacobian) {
            .Object@jacobian.initial <- deriv2(initial, par)
            .Object@jacobian.state <- deriv2(grad, state)
            .Object@jacobian.par <- deriv2(grad, par)
        } else {
            .Object@jacobian.initial <- list()
            .Object@jacobian.state <- list()
            .Object@jacobian.par <- list()
        }

        .Object
    }
)

##' S4 generic for computing a gradient
##' @param object an \code{R} object
##' @param ... further arguments passed to methods
setGeneric(
    "grad",
    def = function(object, ...) {
        standardGeneric("grad")
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
    definition <- function(object, state, par) {
        frame <- as.list(c(state, par))
        gr <- sapply(object@grad, function(grad) eval(grad, as.list(frame)))
        gr
    }
)

##' S4 generic for computing a jacobian
##' @param object an \code{R} object
##' @param ... further arguments passed to methods
setGeneric(
    "jacobian",
    def = function(object, ...) {
        standardGeneric("jacobian")
    }
)

##' Evaluate the jacobian of the gradients
##' @param object object to be evaluated
##' @param state state
##' @param par additional parameters
##' @param type state of par?
##' @docType methods
##' @exportMethod grad
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

setMethod("show", "model.ode",
    function(object){
        cat("Name:", object@name, "\n\n")

        cat("Model:\n")
        f <- paste0("d",object@state, " = ", sapply(object@grad, as.character))
        for(i in 1:length(f))
            cat(f[i], "\n")

        cat("\nInitial values:\n")
        g <- paste0(object@state, "(0) = ", sapply(object@initial, as.character))
        for(i in 1:length(g))
            cat(g[i], "\n")

        cat("\nParameters:", object@par, "\n")
    }
)
