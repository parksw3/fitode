##' S4 generic for evaluating an object
##' @param object an \code{R} object
##' @param ... further arguments passed to methods
##' @return The result of evaluating the object.
##' @keywords internal
##' @export
setGeneric(
    "Eval",
    def = function(object, ...) {
        standardGeneric("Eval")
    }
)

##' S4 generic for computing a gradient
##' @param object an \code{R} object
##' @param ... further arguments passed to methods
##' @return The gradient of the object.
##' @keywords internal
##' @export
setGeneric(
    "grad",
    def = function(object, ...) {
        standardGeneric("grad")
    }
)

##' S4 generic for computing a hessian
##' @param object an \code{R} object
##' @param ... further arguments passed to methods
##' @return The hessian of the object.
##' @keywords internal
##' @export
setGeneric(
    "hessian",
    def = function(object, ...) {
        standardGeneric("hessian")
    }
)

##' S4 generic for computing a jacobian
##' @param object an \code{R} object
##' @param ... further arguments passed to methods
##' @return The jacobian of the object.
##' @keywords internal
##' @export
setGeneric(
    "jacobian",
    def = function(object, ...) {
        standardGeneric("jacobian")
    }
)

##' S4 generic for transforming an object
##' @param object an \code{R} object
##' @param ... further arguments passed to methods
##' @return The result of transforming the object.
##' @keywords internal
##' @export
setGeneric(
    "Transform",
    def = function(object, ...) {
        standardGeneric("Transform")
    }
)
