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
    function(.Object, name,
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
    function(object, transforms=NULL, par) {
        # if no transform, return model
        if (length(transforms) == 0)
            return(object)

        allvars <- c(object@par)
        transforms <- trans(transforms, allvars)

        nstate <- length(object@state)

        newmodel <- newinitial <- vector('list', nstate)

        for(i in 1:nstate) {
            fixed <- c(as.symbol("~"), as.symbol(object@state[i]))
            ff <- lapply(list(object@grad[[i]], object@initial[[i]]), function(x){
                f <- c(fixed, subst(x[[1]]))
                f <- as.formula(as.call(f))
            })
            newmodel[[i]] <- ff[[1]]
            newinitial[[i]] <- ff[[2]]

        }
        if (missing(par)) par <- stop("specify the name of the new parameters")

        new("model.ode", object@name, newmodel, newinitial, par)
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
