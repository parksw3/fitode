to.formula <- function(left, right){
    as.formula(as.call(c(as.symbol("~"), as.symbol(left), right)), env = parent.frame())
}

trans <- function(formulae, allvars) {
    # extract vars from an expression
    vars <- function(e) {
        if (is.numeric(e)) return(c())
        if (is.name(e)) return(as.character(e))
        if (!is.call(e))
            stop("unknown class: ", class(e))
        v <- c()
        for (i in 2:length(e)) {
            v <- c(v, vars(e[[i]]))
        }
        v
    }

    l <- list()
    for (f in formulae) {
        if (!is(f, "formula"))
            stop("transforms must be formula: ", as.character(f))
        var <- as.character(f[[2]])
        if (!var %in% allvars) next
        input <- vars(f[[3]])
        l[[var]] <- f[[3]]
    }
    l
}

# substitute the transformation expressions
subst <- function(e, transforms) {
    if (is.numeric(e)) return(e)
    if (is.name(e)) {
        v <- as.character(e)
        expr <- transforms[[v]]
        if (is.null(expr)) return(e)
        return (expr)
    }
    if (!is.call(e)) stop("unknown class: ", class(e))
    l <- list(e[[1]])
    for (i in 2:length(e))
        l <- c(l, subst(e[[i]], transforms))
    as.call(l)
}

linkfun <- function(link=c("log", "logit")) {
    link <- match.arg(link)
    switch(link,
        log={
            list(
                transform=substitute(exp(x)),
                inverse=substitute(log(x))
            )
        },
        logit={
            list(
                transform=substitute(1/(1+exp(-x))),
                inverse=substitute(log(x)-log(1-x))
            )
        }
    )
}

##' transform parameters
##' @param parms numeric vector containing parameters
##' @param transform list of transformations specified as a formula
##' @param inverse list of inverse transformations specified as formula
##' @examples
##' par <- c(a=2, b=1)
##' transform <- list(a~exp(log.a)+1)
##' inverse <- list(log.a~log(a-1))
##'
##' print(newpar <- transpar(par, transform, inverse))
##'
##' ## direction of transformation is automatically determined
##' print(oldpar <- transpar(newpar, transform, inverse))
##' all.equal(par, oldpar)
##'
##' print(oldpar2 <- transpar(newpar, inverse, transform))
##' all.equal(par, oldpar2)
transpar <- function(parms, transform, inverse) {
    before <- sapply(transform, function(x) as.character(x[[2]]))
    after <- sapply(inverse, function(x) as.character(x[[2]]))

    if (all(before %in% names(parms))) {
        newname <- names(parms)
        newname[match(before, newname)] <- after

        parlist <- lapply(newname, function(x) as.name(x))
        newpar <- lapply(lapply(parlist, subst, trans(inverse,newname)), subst, as.list(parms))
        newpar <- sapply(newpar, eval)
        names(newpar) <- newname
        newpar
    } else if (all(after %in% names(parms))) {
        transpar(parms, inverse, transform)
    } else {
        stop("Wrong transformation/inverse provided?")
    }
}
