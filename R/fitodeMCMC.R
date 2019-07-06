propfun <- function(chol) {
    c(rnorm(ncol(chol), mean=0, sd=1) %*% chol)
}

##' Fit ordinary differential equations model using MCMC
##' @rdname fitodeMCMC
##' @name fitodeMCMC
##' @param model ode model
##' @param data data frame with time column and observation column
##' @param start named vector of starting parameter values
##' @param tcol time column
##' @param vcov
##' @param prior
##' @param mcmc
##' @param burnin
##' @param thin
##' @param link named vector or list of link functions for model parameters
##' @param fixed named vector or list of model parameters to fix
##' @param solver.opts options for ode integration. See \code{\link{ode}}
##' @param solver ode solver
##' @param debug print debugging output?
##' @importFrom mvtnorm dmvnorm
##' @export fitodeMCMC
fitodeMCMC <- function(model, data,
                       start, tcol="times",
                       vcov,
                       prior,
                       mcmc=2000, burnin=1000, thin=1,
                       link,
                       fixed=list(),
                       control=list(maxit=1e5),
                       solver.opts=list(method="rk4"),
                       solver=ode,
                       skip.hessian=FALSE,
                       force.hessian=FALSE,
                       use.ginv=TRUE,
                       ...) {
    if (missing(start)) stop("starting parameters must be specified via 'start'")

    if (missing(vcov)) stop("variance covariance matrix of the proposal distribution must be specified via 'vcov'")

    if (missing(prior)) {
        warning("prior distributions must be specified via 'prior'")
        prior <- function() {return(0)}
    }

    if (length(fixed) > 0) {
        model <- fixpar(model, fixed)
    }

    modelpar <- model@par

    ## check prior
    if (class(prior)=="list") {

    } else if (class(prior)=="function") {

    } else {
        stop("'prior' must be a list of formulas or a function returning log-prior density")
    }



}
