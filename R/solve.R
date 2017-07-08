##' @export
solve <- function(y, times, model, parms,
                 method="lsoda",
                 sensitivity=TRUE,
                 ...) {
    if (sensitivity) {
        nstate <- length(model@state)
        nsens <- nstate * length(model@par)

        ## TODO: allow for expressions in y so that it can be affected by the parameters
        yini <- c(y, rep(0, nsens))

        gfun <- function(times, y, parms) {
            state <- y[1:nstate]
            gr <- grad(model, state, parms)
            js <- jacobian(model, state, parms, type="state")
            jp <- jacobian(model, state, parms, type="par")

            list(c(gr, matrix(y[-c(1:nstate)], ncol=nstate) %*% js + jp))
        }

    } else {
        yini <- y
        gfun <- function(y, times, parms) {
            gr <- grad(model, y, parms)
        }
    }

    ode(yini, times, gfun, parms, method)
}
