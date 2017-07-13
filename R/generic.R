##' S4 generic for evaluate an object
##' @param object an \code{R} object
##' @param ... further arguments passed to methods
setGeneric(
    "Eval",
    def = function(object, ...) {
        standardGeneric("Eval")
    }
)

##' S4 generic for computing a gradient
##' @param object an \code{R} object
##' @param ... further arguments passed to methods
setGeneric(
    "grad",
    def = function(object, ...) {
        standardGeneric("grad")
    }
)

##' S4 generic for computing a hessian
##' @param object an \code{R} object
##' @param ... further arguments passed to methods
setGeneric(
    "hessian",
    def = function(object, ...) {
        standardGeneric("hessian")
    }
)

##' S4 generic for computing a jacobian
##' @param object an \code{R} object
##' @param ... further arguments passed to methods
setGeneric(
    "jacobian",
    def = function(object, ...) {
        standardGeneric("jacobian")
    }
)

##' S4 generic for transforming an object
##' @param object an \code{R} object
##' @param ... further arguments passed to methods
setGeneric(
    "Transform",
    def = function(object, ...) {
        standardGeneric("Transform")
    }
)