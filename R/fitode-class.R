##' Class "odemodel" representing ode models
##'
##' @slot name name of the model
##' @slot gfun gradient function
##' @slot grad list of gradients
##' @slot observation list of observation models
##' @slot initial list of expressions representing the initial values
##' @slot state state variables
##' @slot par parameters
##' @slot link link functions for parameters (log links are used as default)
##' @slot diffnames character vector specifying the names of a variable for which the consecutive difference needs to be calculated
##' @slot call original function call
##' @slot jacobian.initial Jacobian of initial values with respect to its parameters
##' @slot jacobian.state Jacobian with respect to its states
##' @slot jacobian.par Jacobian with repsect to its parameters
##' @slot loglik list of log-likelihood functions
##' @slot expr expressions for true trajectories
##' @slot expr.sensitivity sensitivity of the expressions with respect to state variables and parameters
##' @slot keep_sensitivity (logical) keep sensitivity equations?
##' @exportClass odemodel
setClass(
    "odemodel",
    slots = c(
        name = "character",
        gfun = "function",
        grad = "list",
        observation = "list",
        initial= "list",
        state = "character",
        par = "character",
        link = "character",
        diffnames = "character",
        jacobian.initial = "list",
        jacobian.state = "list",
        jacobian.par = "list",
        loglik = "list",
        expr="list",
        expr.sensitivity="list",
        keep_sensitivity  = "logical",
        call = "language"
    )
)

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

##' Class representing prior models used to fit ode models
##'
##' @slot name name of the distribution
##' @slot expr an expression specifying the model
##' @slot observation observation variable name
##' @slot par additional parameter names
##' @slot keep_grad keep gradient?
##' @slot grad the gradient with respect to the parameters
##' @exportClass prior.ode
setClass(
    "prior.ode",
    slots = c(
        name = "character",
        expr = "expression",
        observation = "character",
        par = "character",
        keep_grad = "logical",
        grad = "list"
    )
)

##' Class "fitode".
##' Result of ode fitting based on Maximum Likelihood Estimation
##'
##' @slot call (languge) The call to \code{\link{fitode}}
##' @slot model odemodel object
##' @slot data the time series data
##' @slot coef estimated parameters
##' @slot vcov estimated variance-covariance matrix
##' @slot min minimum negative log-likelihood
##' @slot mle2 mle2 object
##' @slot link list of link functions for model parameters
##' @slot fixed list of fixed parameters
##' @slot prior list of priors
##' @seealso \code{\link{mle2-class}}
##' @exportClass fitode
setClass("fitode",
    slots = c(
        call="language",
        model="odemodel",
        data="data.frame",
        coef="numeric",
        vcov="matrix",
        min="numeric",
        mle2="mle2",
        link="list",
        fixed="list",
        prior="list"
    )
)

## need this to use mcmc.list in S4?
setOldClass("mcmc.list")

##' Class "fitodeMCMC".
##' Result of ode fitting based on Markov Chain Monte Carlo (MCMC)
##'
##' @slot call (languge) The call to \code{\link{fitodeMCMC}}
##' @slot model odemodel object
##' @slot data the time series data
##' @slot coef estimated parameters (posterior median)
##' @slot vcov estimated variance-covariance matrix
##' @slot mcmc mcmc.list object containing posterior samples
##' @slot lp mcmc.list object containing log-posterior values of posterior samples
##' @slot link list of link functions for model parameters
##' @slot fixed list of fixed parameters
##' @slot prior list of priors
##' @slot details a list containing miscellaneous objects for internal uses
##' @exportClass fitodeMCMC
setClass("fitodeMCMC",
         slots = c(
             call="language",
             model="odemodel",
             data="data.frame",
             coef="numeric",
             vcov="matrix",
             mcmc="mcmc.list",
             lp="mcmc.list",
             link="list",
             fixed="list",
             prior="list",
             details="list"
         )
)
