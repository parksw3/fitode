## get LHS of formula (element 3), then take its head (element 1)
get_head <- function(x) deparse(get_rhs(x)[[1]])

get_rhs <- function(x) x[[3]]

## copied from lme4::namedList
named_list <- function (...) {
    L <- list(...)
    snm <- sapply(substitute(list(...)), deparse)[-1]
    if (is.null(nm <- names(L))) 
        nm <- snm
    if (any(nonames <- nm == "")) 
        nm[nonames] <- snm[nonames]
    setNames(L, nm)
}
