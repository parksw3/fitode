## get LHS of formula (element 3), then take its head (element 1)
get_head <- function(x) deparse(get_rhs(x)[[1]])

get_rhs <- function(x) x[[3]]
