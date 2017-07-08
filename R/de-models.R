##' Class representing ode models
##'
##' @slot name name of the model
##' @slot model list of formulas specifying ode model
##' @slot gradient ode gradients
##' @slot state state variables
##' @slot par parameters
##' @slot js Jacobian with respect to its states
##' @slot jp Jacobian with repsect to its parameters
setClass(
    "de-model",
    slots = c(
        name = "character",
        model = "list",
        state = "character",
        par = "character",
        gradient = "list",
        js = "list",
        jp = "list"
    )
)

##' the initializer for de-model
##'
##' @param .Object object
##' @slot name name of the model
##' @slot model ode model
##' @slot state state variables
##' @slot par parameters
##' @slot keep_jacobian (logical) maintain the Jacobian as part of the model
##' @examples
##' SI_model <- new("de-model",
##'     name = "SI",
##'     model = list(
##'         S ~ - beta*S*I,
##'         I ~ beta*S*I - gamma*I
##'     ),
##'     state = c("S", "I"),
##'     par= c("beta", "gamma")
##' )
##' @docType methods
##' @exportMethod initialize
setMethod(
    "initialize",
    "de-model",
    definition = function(.Object, name,
                          model,
                          state, par,
                          keep_jacobian=TRUE) {
        .Object@name <- name
        if (any(unlist(lapply(model, class)) != "formula"))
            stop("model must be a list of formulas")

        if (length(model) != length(state))
            stop("model must have the same length as the number of states")

        ## a bit awkward
        if (!all.equal(as.character(lapply(model, "[[", 2)), state))
            stop("model does not match the states provided?")

        .Object@model <- model
        .Object@state <- state
        .Object@par <- par

        gradient <- lapply(model, function(x) as.expression(x[[3]]))

        .Object@gradient <- gradient

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
            names(d) <- vars
            d
        }

        if (keep_jacobian) {
            .Object@js <- deriv2(gradient, state)
            .Object@jp <- deriv2(gradient, par)
        } else {
            .Object@js <- list()
            .Object@jp <- list()
        }

        .Object
    }
)

setMethod("show", "de-model",
    function(object){
        cat("Name:", object@name, "\n\n")
        lapply(object@model, function(x){
            f <- paste(x[2], x[3], sep=" = ")
            f <- paste0("d", f)
            cat(f, "\n")
        })
        cat("\nStates:", object@state, "\n")
        cat("\nParameters:", object@par, "\n")
    }
)


