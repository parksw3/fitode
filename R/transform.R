#' @importFrom methods is
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
