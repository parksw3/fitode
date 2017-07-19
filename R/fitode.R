##' fit ode
##' @param formula formula specifing observation variable and mean
##' @param data data frame
##' @import bbmle
##' @export
fitode <- function(formula, start,
                   model, loglik,
                   data,
                   tcol = "times",
                   control=list(maxit=1e5),
                   debug=FALSE) {
    ocol <- as.character(formula[[2]])
    data <- data.frame(times = data[[tcol]], observation = data[[ocol]])

    dataarg <- c(data,list(model=model, loglik=loglik, formula=formula))

    f.env <- new.env()
    ## set initial values
    assign("oldnll",NULL,f.env)
    assign("oldpar",NULL,f.env)
    assign("oldgrad",NULL,f.env)
    assign("data", data, f.env)

    objfun <- function(par, formula, model, loglik, observation, times) {
        if (identical(par,oldpar)) {
            if (debug) cat("returning old version of value\n")
            return(oldnll)
        }
        if (debug) cat("computing new version (nll)\n")

        v <- ode.sensitivity(par, formula, model, loglik, observation, times)
        oldnll <<- v[1]
        oldgrad <<- v[-1]
        oldpar <<- par

        return(oldnll)

    }
    gradfun <- function(par, formula, model, loglik, observation, times) {
        if (identical(par,oldpar)) {
            if (debug) cat("returning old version of grad\n")
            return(oldgrad)
        }
        if (debug) cat("computing new version (grad)\n")

        v <- ode.sensitivity(par, formula, model, loglik, observation, times)
        oldnll <<- v[1]
        oldgrad <<- v[-1]
        oldpar <<- par
        return(oldgrad)
    }

    environment(objfun) <- f.env

    environment(gradfun) <- f.env

    parnames <- c(model@par, loglik@par)
    attr(objfun, "parnames") <- parnames

    m <- mle2(objfun,
              vecpar=TRUE,
              start=start,
              method="BFGS",
              control=control,
              gr=gradfun,
              data=dataarg)
    return(m)
}

ode.sensitivity <- function(parms, formula,
                        model, loglik,
                        observation, times=NULL) {
    if (is.null(times)) times <- seq(length(count))
    solution <- solve(model, times, parms)
    expr <- as.expression(formula[[3]])

    mean <- eval(expr, solution@solution) ## TODO: what happens if mean is affected by the parameters?

    nstate <- length(model@state)

    dm <- lapply(model@state, function(s) Deriv(expr, s))
    sens <- vector('list', nstate)
    for(i in 1:nstate) {
        sens[[i]] <- eval(dm[[i]], solution@solution) * solution@sensitivity[[i]]
    }

    sens <- do.call("+", sens)

    loglik.par <- as.list(parms[-c(1:length(model@par))])

    nll <- -sum(Eval(loglik, observation, mean, loglik.par))
    loglik.gr <- grad(loglik, observation, mean, loglik.par)
    sensitivity <- c(-colSums(loglik.gr[[1]] * sens))
    if(length(loglik.gr) > 1) sensitivity <- c(sensitivity, -sapply(loglik.gr[-1], sum))

    c(nll, sensitivity)
}
