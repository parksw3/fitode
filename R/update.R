##' Update fitode fits
##' @param object fitode objects
##' @param observation observation model
##' @param initial initial values
##' @param par model parameters
##' @param link link functions for parameters (log links are used as default)
##' @param ... additional arguments to be passed to fitode
##' @exportMethod update
setMethod("update", "fitode",
    function (object,
              observation, initial,
              par, link,
              ...){
        call <- object@call

        extras <- match.call(expand.dots = FALSE)$...

        model <- eval(call$model)

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

        eval(call, parent.frame())
    }
)
