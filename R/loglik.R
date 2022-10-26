##' The initializer for loglik.ode
##'
##' @param .Object object
##' @param name name of the distribution
##' @param model formula specifying the disstribution
##' @param observation observation variable name
##' @param mean mean variable name
##' @param par additional parameter names
##' @param keep_grad maintain the gradient as part of the model
##' @return An object of class ``loglik.ode'' as described in \code{\link{loglik.ode-class}}.
##' @importFrom Deriv Deriv
##' @docType methods
##' @keywords internal
##' @exportMethod initialize
setMethod(
    "initialize",
    "loglik.ode",
    function(.Object, name, model, observation="X", mean, par=NULL,
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

##' Evaluate the log-likelihood model
##' @param object loglik.ode object
##' @param observation observations
##' @param mean mean values
##' @param par additional parameters
##' @param ... other values if required
##' @return numeric
##' @docType methods
##' @keywords internal
##' @exportMethod Eval
setMethod(
    "Eval",
    "loglik.ode",
    function(object, observation, mean, par=NULL, ...) {
        frame <- list(observation, mean)
        frame <- append(frame, par)
        names(frame) <- c(object@observation, object@mean, object@par)
        frame <- append(frame, list(...))
        eval(object@expr, frame)
    }
)

##' Evaluate the gradient of a likelihood model
##' @param object loglik.ode object
##' @param observation observations
##' @param mean mean values
##' @param par additional parameters
##' @param ... other values if required
##' @return a list with each element as a partial derivative values
##' @docType methods
##' @keywords internal
##' @exportMethod grad
setMethod(
    "grad",
    "loglik.ode",
    function(object, observation, mean, par, ...) {
        frame <- list(observation, mean)
        frame <- append(frame, par)
        names(frame) <- c(object@observation, object@mean, object@par)
        frame <- append(frame, list(...))
        l <- lapply(object@grad, function(deriv) { eval(deriv, frame)})
        l
    }
)

##' Transform the model
##' @param object loglik.ode object
##' @param name name of the log-likelihood model
##' @param transforms list of formulas specifying transformations
##' @param observation observation variable name
##' @param mean mean variable name
##' @param par additional parameter names
##' @param keep_grad maintain the gradient as part of the model
##' @return loglik.ode object
##' @docType methods
##' @keywords internal
##' @exportMethod Transform
setMethod(
    "Transform",
    "loglik.ode",
    function(object, transforms=NULL,
             name,
             observation="X",
             mean, par,
             keep_grad=TRUE) {
        # if no transform, return model
        if (length(transforms) == 0)
            return(object)
        allvars <- c(object@observation, object@mean, object@par)
        transforms <- trans(transforms, allvars)
        f <- c(as.symbol("~"), as.symbol("LL"), subst(object@expr[[1]], transforms))
        f <- as.formula(as.call(f))

        if (missing(name)) name <- object@name
        if (missing(mean)) mean <- object@mean
        if (missing(par)) par <- object@par

        new("loglik.ode", name, f, observation, mean=mean, par=par, keep_grad=keep_grad)
    }
)


## discontinuity in second derivative, but ... probably OK
## TODO: not sure why it won't work unless I export it??
##' Taylor expansion of digamma(a+b) for a>>b
##' @param x first argument
##' @param y second argument
##' @param mag cutoff magnitude for switching approximations
##' @return numeric
##' @importFrom Deriv drule
##' @keywords internal
##' @export
dfun <- function(x,y,mag=1e8) {
    return(ifelse(x/y>mag,
                  -y*trigamma(x),
                  digamma(x)-digamma(x+y)))
}

##' Taylor expansion of trigamma(a+b) (?) for a>>b
##' @param x first argument
##' @param y second argument
##' @param mag cutoff magnitude for switching approximations
##' @return numeric
##' @keywords internal
##' @export
dfun2 <- function(x,y,mag=1e8,focal="x") {
    return(switch(focal,
                  x=ifelse(x/y>mag,
                           -y*psigamma(x,2),
                           trigamma(x)-trigamma(x+y)),
                  y=ifelse(x/y>mag,
                           -trigamma(x),
                           -trigamma(x+y))))
}

NBconst <- function(k,x) {
    return(ifelse(x==0,0,lbeta(k,x)+log(x)))
}

## TODO: call this family instead of dist... but probably doesn't matter
## for internal usage

##' Select a log-likelihood model
##' @keywords internal
##' @param dist conditional distribution of reported data (dnorm, dnorm2, dpois, dnbinom, dnbinom1, dgamma)
select_model <- function(dist = c("ols", "dnorm", "dnorm2", "dpois", "dnbinom", "dnbinom1", "dgamma", "dlnorm")) {
    dist <- match.arg(dist)
    name <- dist
    model <- switch(dist,
        ## FIXME: in general, why create a named object then return it by name ... ?
        ols={
            loglik_ols <- new("loglik.ode", "ols",
                LL ~ -(X-mean)^2,
                mean="mean",par=NULL)
            loglik_ols
        }, dnorm={
            loglik_gaussian <- new("loglik.ode", "gaussian",
                LL ~ -(X-mean)^2/(2*sd^2) - log(sd) - 1/2*log(2*pi),
                mean="mean", par="sd")

            loglik_gaussian
        }, dnorm2={
            loglik_gaussian2 <- new("loglik.ode", "gaussian",
                                   LL ~ -(X-mean)^2/(2*sd^2) - log(sd) - 1/2*log(2*pi),
                                   mean="mean", par="sd")

            loglik_gaussian2 <- Transform(loglik_gaussian2,
                                         transforms = list(sd ~ sqrt(sum((X-mean)^2)/(length(X)-1))),
                                         par=NULL)

            loglik_gaussian2
        }, dpois={
            loglik_poisson <- new("loglik.ode", "poisson",
                LL ~ X*log(lambda) - lambda - lgamma(X+1),
                mean = "lambda", par = c())
            loglik_poisson
        }, dnbinom={
            loglik_nbinom <- new ("loglik.ode", "nbinom",
                LL ~ - NBconst(size, X) + size * (-log1p(mu/size)) +
                    X * log(mu) - X * log(size + mu),
                mean="mu",
                par = "size")

            loglik_nbinom
        }, dnbinom1={
            loglik_nbinom1 <- new ("loglik.ode", "nbinom1",
                LL ~ - NBconst(mu/phi, X) + mu/phi * (-log1p(phi)) +
                    X * log(mu) - X * log(mu/phi + mu),
                mean="mu",
                par = "phi")

            loglik_nbinom1
        }, dgamma={
            loglik_gamma <- new("loglik.ode", "gamma",
                LL ~ shape * log(shape/mean) + (shape - 1) * log(X) - shape/mean * X - lgamma(shape),
                mean="mean",
                par="shape"
            )

            loglik_gamma
        },
        dlnorm={
            loglik_lognormal <- new("loglik.ode", "lognormal",
                LL ~-(log(X)-meanlog)^2/(2*sdlog^2) -
                    log(sdlog) - 1/2*(log(2*pi)) -log(sdlog) - log(X),
                ## 1/(sqrt(2 pi) sigma x) e^-((log x - mu)^2 / (2 sigma^2))
                mean="meanlog", par="sdlog")
            loglik_lognormal
    })

    model@name <- name

    model
}
