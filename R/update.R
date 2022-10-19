##' Update fitode fits
##' @param object fitode objects
##' @param observation observation model
##' @param initial initial values
##' @param par model parameters
##' @param link link functions for parameters (log links are used as default)
##' @param ... additional arguments to be passed to fitode
##' @return An object of class ``fitode'' as described in \code{\link{fitode-class}}.
##' @exportMethod update
setMethod("update", "fitode",
    function (object,
              observation, initial,
              par, link,
              ...){
        pf <- parent.frame(2) ## need this to jump back through the S4 '.local' context?
        ## pass this down to avoid too much guessing about how far back up to go ...
        update_internal(object, observation, initial, par, link, pf= pf, ...)
    }
)

##' Update fitodeMCMC fits
##' @param object fitodeMCMC objects
##' @param observation observation model
##' @param initial initial values
##' @param par model parameters
##' @param link link functions for parameters (log links are used as default)
##' @param ... additional arguments to be passed to fitode
##' @return An object of class ``fitode'' as described in \code{\link{fitodeMCMC-class}}.
##' @exportMethod update
setMethod("update", "fitodeMCMC",
    function (object,
              observation, initial,
              par, link,
              ...){
        update_internal(object, observation, initial, par, link, pf=parent.frame(), ...)
    }
)

update_internal <- function(object,
                            observation, initial,
                            par, link,
                            pf,
                            ...){
    call <- object@call
    ## FIXME: why doesn't think work?
    ## extras <- match.call(expand.dots = FALSE)$...
    extras <- list(...)
    model <- eval(call$model, pf)

    if (!missing(observation) || !missing(initial) || !missing(par) || !missing(link)) {
        model <- Transform(
            model,
            observation=observation,
            initial=initial,
            par=par,
            link=link
        )

        call$model <- model
    }

    ## taken from bbmle
    ## https://github.com/bbolker/bbmle/blob/master/R/update.R
    if (length(extras)) {
        existing <- !is.na(match(names(extras), names(call)))
        for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
        if (any(!existing)) {
            call <- c(as.list(call), extras[!existing])
            call <- as.call(call)
        }
    }

    ## CHECK:
    ## I probably have to go up twice to evaluate this properly?
    ## bbmle::update() goes up once but fitode has update_internal
    eval(call, pf)
}
