##' Class representing ode models
##'
##' @slot name name of the model
##' @slot gfun gradient function
##' @slot initial list of expressions representing the initial values
##' @slot state state variables
##' @slot par parameters
##' @slot jacobian.initial Jacobian of initial values with respect to its parameters
##' @slot jacobian.state Jacobian with respect to its states
##' @slot jacobian.par Jacobian with repsect to its parameters
##' @slot keep_sensitivity (logical) keep sensitivity equations?
##' @exportClass model.ode
setClass(
    "model.ode",
    slots = c(
        name = "character",
        gfun = "function",
        initial= "list",
        state = "character",
        par = "character",
        jacobian.initial = "list",
        jacobian.state = "list",
        jacobian.par = "list",
        keep_sensitivity  = "logical"
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

##' Class "fitode".
##' Result of ode fitting based on Maximum Likelihood Estimation
##'
##' @name fitode
##' @rdname fitode
##' @seealso \code{\link{mle2-class}}
##' @exportClass fitode
##' @export
fitode <- setClass("fitode",
    slots = c(
        formula="formula",
        model="model.ode",
        loglik="loglik.ode",
        data="data.frame",
        coef="numeric",
        vcov="matrix",
        min="numeric",
        mle2="mle2",
        link="list"
    )
)


