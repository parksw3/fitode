##' the initializer for prior.ode
##'
##' @param .Object object
##' @param name name of the distribution
##' @param model the formula specifying the model
##' @param observation observation variable name
##' @param par additional parameter names
##' @param keep_grad maintain the gradient as part of the model
##' @return An object of class ``prior.ode'' as described in \code{\link{prior.ode-class}}.
##' @docType methods
##' @keywords internal
##' @exportMethod initialize
setMethod(
    "initialize",
    "prior.ode",
    function(.Object, name, model, observation="X", par=NULL,
             keep_grad=TRUE) {
        .Object@name <- name
        if (!is(model, "formula"))
            stop("model must be a formula")
        f <- as.list(model)
        .Object@expr <- as.expression(f[[3]])
        .Object@observation <- observation
        .Object@par <- par <- as.character(par)

        if(.Object@keep_grad <- keep_grad) {
            # compute the gradient
            # with respect to X this time
            vars <- observation
            deriv <- function(expr) {
                d <- lapply(vars,
                            function(p){
                                ## hack to make drule available in current environment
                                Deriv(expr, p)
                            })
                names(d) <- vars
                d
            }

            .Object@grad <- deriv(.Object@expr)
        } else .Object@grad <- list()

        .Object
    }
)

##' Transform the prior model
##' @param object object
##' @param name name of the log-likelihood model
##' @param transforms list of formulas specifying transformations
##' @param observation observation variable name
##' @param par additional parameter names
##' @param keep_grad maintain the gradient as part of the model
##' @return An object of class ``prior.ode'' as described in \code{\link{prior.ode-class}}.
##' @docType methods
##' @keywords internal
##' @exportMethod Transform
setMethod(
    "Transform",
    "prior.ode",
    function(object, transforms=NULL,
             name,
             observation="X",
             par,
             keep_grad=TRUE) {
        # if no transform, return model
        if (length(transforms) == 0)
            return(object)
        allvars <- c(object@observation, object@par)
        transforms <- trans(transforms, allvars)
        f <- c(as.symbol("~"), as.symbol("LL"), subst(object@expr[[1]], transforms))
        f <- as.formula(as.call(f))

        if (missing(name)) name <- object@name
        if (missing(par)) par <- object@par
        if (missing(keep_grad)) keep_grad <- object@keep_grad

        new("prior.ode", name, f, observation, par=par, keep_grad=keep_grad)
    }
)

##' Make a list containing log prior density and its gradient
##' @param model model.ode object
##' @param link link
##' @param prior list of formulas
##' @param prior.density (logical) does this represent a probability density?
##' @param keep_grad (logical) keep gradient?
##' @keywords internal
make_prior <- function(model, link, prior, prior.density=TRUE, keep_grad=TRUE) {
    prior_list <- lapply(prior, function(ll) {
        tlink <- link[match(deparse(ll[[2]]), names(link))]

        ll_model <- select_prior(as.character(ll[[3]][[1]]), link=tlink, prior.density=prior.density, keep_grad=keep_grad)

        tpar <- ifelse(tlink=="identity", deparse(ll[[2]]), paste0(tlink, ".", deparse(ll[[2]])))

        trans_obs <- as.formula(as.call(c(as.symbol("~"), as.symbol("X"), as.symbol(tpar))))

        trans_list <- list(trans_obs)

        likpar <- ll_model@par

        if(!all(likpar %in% names(as.list(ll[[3]]))))
            stop(paste0("'", deparse(ll[[3]][[1]]), "' must be parameterized with ", paste(paste0("'", likpar, "'"), collapse=" and ")))

        call <- as.list(ll[[3]])[likpar]

        mm <- Map(function(x, y) as.formula(as.call(c(as.symbol("~"), as.symbol(x), y))), x=likpar, y=call)

        trans_list <- append(trans_list, mm)

        ll_model <- Transform(ll_model,
                              observation=tpar,
                              transforms=trans_list,
                              par=character(0))

        ll_model
    })

    gradlist <- vector('list', length(model@par))
    names(gradlist) <- model@par

    if (keep_grad)
        gradlist[sapply(prior_list, function(x) names(x@observation))] <- sapply(prior_list, function(x) x@grad)

    list(
        prior.density=parse(text=paste(sapply(prior_list, function(x) as.character(x@expr)), collapse="+")),
        prior.grad=gradlist
    )
}

## very similar to select_model

##' Select a prior model
##' @param family prior distribution type
##' @param link link
##' @param prior.density (logical) keep the Jacobian of transformations?
##' @param keep_grad (logical) keep gradients?
##' @keywords internal
select_prior <- function(family = c("dnorm", "dgamma", "dbeta",
                                    "dlnorm"),
                         link = c("identity", "log", "logit"),
                         prior.density=TRUE,
                         keep_grad=TRUE) {
    family <- match.arg(family)
    link <- match.arg(link)
    name <- family

    if (family=="dbeta" && link=="log")
        stop("beta prior distribution is incompatible with log link")

    if (family!="dbeta" && link=="logit")
        stop("logit link can be only used with beta prior distribution")


    model <- switch(family,
        dnorm={
            loglik_gaussian <- new("prior.ode", "gaussian",
                                   LL ~ -(X-mean)^2/(2*sd^2) - log(sd) - 1/2*log(2*pi) + constant,
                                   par=c("mean", "sd", "constant"),
                                   ## "constant" is the d(mu)/d(eta) transformation
                                   keep_grad=keep_grad)
            loglik_gaussian
        }, dgamma={
            loglik_gamma <- new("prior.ode", "gamma",
                                LL ~ shape * log(rate) + (shape - 1) * log(X) - rate * X - lgamma(shape) + constant,
                                par=c("shape", "rate", "constant"),
                                keep_grad=keep_grad)

            loglik_gamma
        }, dbeta={
            loglik_beta <- new("prior.ode", "beta",
                               LL ~ (shape1-1) * log(X) + (shape2-1) * log(1-X) -
                                   lgamma(shape1) - lgamma(shape2) + lgamma(shape1 + shape2) + constant,
                               par=c("shape1", "shape2", "constant"),
                               keep_grad=keep_grad)
        },
        dlnorm={
            loglik_lognormal <- new("prior.ode", "lognormal",
                           LL ~-(log(X)-meanlog)^2/(2*sdlog^2) -
                              log(sdlog) - 1/2*(log(2*pi)) -log(sdlog) - log(X),
                ## 1/(sqrt(2 pi) sigma x) e^-((log x - mu)^2 / (2 sigma^2))
                par = c("meanlog", "sdlog", "constant"),
                keep_grad=keep_grad)
            loglik_lognormal
        }
    )

    model@name <- family

    if (link=="log") {
        model <- Transform(
            model,
            transforms=list(X~exp(X))
        )
    } else if (link=="logit") {
        model <- Transform(
            model,
            transforms=list(X~exp(X)/(1 + exp(X)))
        )
    }

    if (prior.density) {
        if (link=="log") {
            model <- Transform(
                model,
                transforms=list(constant ~ X),
                par=head(model@par, -1)
            )
        } else if (link=="logit") {
            model <- Transform(
                model,
                transforms=list(constant ~ X - 2 * log(exp(X) + 1)),
                par=head(model@par, -1)
            )
        }
    } else {
        model <- Transform(
            model,
            transforms=list(constant ~ 0),
            par=head(model@par, -1)
        )
    }

    model
}
