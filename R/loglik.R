##' Class representing log-likelihood models used to fit ode models
##'
##' @slot name name of the distribution
##' @slot expr an expression specifying the model
##' @slot observation observation variable name
##' @slot mean mean variable name
##' @slot par additional parameter names
##' @slot grad the gradient with respect to the parameters
##' @exportClass loglik.ode
setClass(
    "loglik.ode",
    slots = c(
        name = "character",
        expr = "expression",
        observation = "character",
        mean = "character",
        par = "character",
        grad = "list"
    )
)

##' the initializer for loglik.ode
##'
##' @param .Object object
##' @param name name of the distribution
##' @param model the formula specifying the model
##' @param observation observation variable name
##' @param mean mean variable name
##' @param par additional parameter names
##' @param keep_grad maintain the gradient as part of the model
##' @importFrom Deriv Deriv
##' @docType methods
##' @exportMethod initialize
setMethod(
    "initialize",
    "loglik.ode",
    definition = function(.Object, name, model, observation="X", mean, par=NULL,
                          keep_grad=TRUE) {
        .Object@name <- name
        if (!is(model, "formula"))
            stop("model must be a formula")
        f <- as.list(model)
        .Object@expr <- as.expression(f[[3]])
        .Object@observation <- observation
        .Object@mean <- mean
        .Object@par <- par <- as.character(par)
        # compute the gradient
        vars <- c(mean, par)
        deriv <- function(expr) {
            d <- lapply(vars,
                        function(p){
                            ## hack to make drule available in current environment
                            Deriv(expr, p)
                        })
            names(d) <- vars
            d
        }
        if(keep_grad) {
            .Object@grad <- deriv(.Object@expr)
        } else .Object@grad <- list()

        .Object
    }
)

##' Evaluate log likelihood model
##' @param object object to be evaluated
##' @param observation observations
##' @param mean mean values
##' @param par additional parameters
##' @param ... other values if required
##' @return numeric
##' @docType methods
##' @exportMethod Eval
setMethod(
    "Eval",
    "loglik.ode",
    definition = function(object, observation, mean, par=NULL, ...) {
        frame <- list(observation, mean)
        frame <- append(frame, par)
        names(frame) <- c(object@observation, object@mean, object@par)
        frame <- append(frame, list(...))
        eval(object@expr, frame)
    }
)

##' Evaluate the gradient of a likelihood model
##' @param object object to be evaluated
##' @param observation observations
##' @param mean mean values
##' @param par additional parameters
##' @param ... other values if required
##' @return a list with each element as a partial derivative values
##' @docType methods
##' @exportMethod grad
setMethod(
    "grad",
    "loglik.ode",
    definition <- function(object, observation, mean, par, ...) {
        frame <- list(observation, mean)
        frame <- append(frame, par)
        names(frame) <- c(object@observation, object@mean, object@par)
        frame <- append(frame, list(...))
        l <- lapply(object@grad, function(deriv) { eval(deriv, frame)})
        l
    }
)

##' Transform the model
##' @param object object
##' @param name name of the log-likelihood model
##' @param transforms list of formulas specifying transformations
##' @param observation observation variable name
##' @param mean mean variable name
##' @param par additional parameter names
##' @param keep_grad maintain the gradient as part of the model
##' @return loglik.ode object
##' @docType methods
##' @exportMethod Transform
setMethod(
    "Transform",
    "loglik.ode",
    definition <- function(object, transforms=NULL,
                           name,
                           observation="X",
                           mean, par,
                           keep_grad=TRUE) {
        # if no transform, return model
        if (length(transforms) == 0)
            return(object)
        allvars <- c(object@observation, object@mean, object@par)
        transforms <- trans(transforms, allvars)
        f <- to.formula("LL", subst(object@expr[[1]], transforms))

        if (missing(name)) name <- object@name
        if (missing(mean)) mean <- object@mean
        if (missing(par)) par <- object@par

        new("loglik.ode", name, f, observation, mean=mean, par=par, keep_grad=keep_grad)
    }
)

## use Taylor expansion of digamma(a+b) for a>>b
## discontinuity in second derivative, but ... probably OK
##' @importFrom Deriv drule
dfun <- function(x,y,mag=1e8) {
    return(ifelse(x/y>mag,
                  -y*trigamma(x),
                  digamma(x)-digamma(x+y)))
}
dfun2 <- function(x,y,mag=1e8,focal="x") {
    return(switch(focal,
                  x=ifelse(x/y>mag,
                           -y*psigamma(x,2),
                           trigamma(x)-trigamma(x+y)),
                  y=ifelse(x/y>mag,
                           -trigamma(x),
                           -trigamma(x+y))))
}

w_lbeta <- function(a,b) {
    ## when we have an effectively-Poisson case
    ## lbeta gives "underflow occurred in 'lgammacor'" frequently ...
    ## suppressWarnings() causes an obscure error ?
    ## using w_lbeta rather than lbeta causes obscure errors from Deriv()
    op <- options(warn=-1)
    on.exit(options(op))
    return(lbeta(a,b))
}

##' Select likelihood model
##' @param dist conditional distribution of reported data
##' @export
select_model <- function(dist = c("gaussian", "poisson", "quasipoisson", "nbinom", "nbinom1")) {
        dist <- match.arg(dist)
    name <- dist
    if (dist == "quasipoisson") dist <- "poisson"
    model <- switch(dist,
        gaussian={
            loglik_gaussian <- new("loglik.ode", "gaussian",
                LL ~ -(X-mu)^2/(2*sigma^2) - log(sigma) - 1/2*log(2*pi),
                mean="mu", par="sigma")

            loglik_gaussian
        }, poisson={
            loglik_poisson <- new("loglik.ode", "poisson",
                LL ~ X*log(lambda) - lambda - lgamma(X+1),
                mean = "lambda", par = c())
            loglik_poisson
        }, nbinom={
            loglik_nbinom <- new ("loglik.ode", "nbinom",
                LL ~ -lbeta(k, X) - log(X) + k * (-log1p(mu/k)) +
                    X * log(mu) - X * log(k + mu),
                mean="mu",
                par = "k")

            loglik_nbinom <- Transform(
                loglik_nbinom,
                transforms = list(k ~ exp(ll.k)),
                par="ll.k"
            )

            loglik_nbinom
        }, nbinom1={
            loglik_nbinom1 <- new ("loglik.ode", "nbinom",
                LL ~ -lbeta(mu/phi, X) - log(X) + mu/phi * (-log1p(phi)) +
                    X * log(mu) - X * log(mu/phi + mu),
                mean="mu",
            par = "phi")

            loglik_nbinom1 <- Transform(
                loglik_nbinom1,
                    transforms = list(phi ~ exp(ll.phi)),
                par=c("ll.phi")
            )

            loglik_nbinom1
        }
    )

    model@name <- name

    model
}
