##' fit ode
##' @rdname fitode
##' @name fitode
##' @param formula formula specifing observation variable and mean
##' @param start starting values for optimization
##' @param model ode model
##' @param loglik log liklihood model
##' @param data data frame with time column and observation column
##' @param tcol time column
##' @param links named vector or list of link functions for ode parameters
##' @param control see optim
##' @param ode.opts options for ode integration. See ode
##' @param debug print debugging output?
##' @import bbmle
##' @importFrom numDeriv jacobian
##' @export fitode
setMethod(
    "initialize",
    "fitode",
    function(.Object,
             formula, start,
             model, loglik=select_model("gaussian"),
             data,
             tcol = "times",
             links,
             control=list(maxit=1e5),
             ode.opts=list(method="rk4", hini=0.1),
             debug=FALSE) {
        if (any(is.na(match(names(start), c(model@par, loglik@par))))) {
            stop(
                paste0("`start` must specify the following parameters:\n",
                    "\node parameters: ", paste(model@par, collapse = ", "),
                    "\nlikelihood parameters: ", paste(loglik@par)
                )
            )
        }


        orig.model <- model
        .Object@model <- orig.model
        orig.formula <- formula
        .Object@formula <- orig.formula
        .Object@loglik <- loglik

        ocol <- as.character(formula[[2]])
        data <- data.frame(times = data[[tcol]], observation = data[[ocol]])

        .Object@data <- data

        if (!missing(links)) {
            .Object@links <- links

            if(!is.list(links)) links <- as.list(links)

            transform <- vector('list', length(links))
            inverse <- vector('list', length(links))

            newpar <- oldpar <- model@par

            for(i in 1:length(links)) {
                par <- names(links)[i]
                link <- links[[i]]

                tpar <- paste(link, par, sep=".")

                newpar[which(newpar==par)] <- tpar

                m <- Map(subst, e=linkfun(link), transforms=list(list(x=as.name(tpar)), list(x=as.name(par))))

                transform[[i]] <- to.formula(par, m$transform)
                inverse[[i]] <- to.formula(tpar, m$inverse)
            }

            .Object@transforms <- list(
                transform=transform,
                inverse=inverse
            )

            model <- Transform(model, transforms=transform, par=newpar)

            newformula <- subst(formula[[3]], trans(transform, oldpar))
            formula <- to.formula(formula[[2]], newformula)

            start <- transpar(start, transform, inverse)
        } else {
            .Object@links <- list()
            .Object@transforms <- list()
        }

        dataarg <- c(data,list(model=model, loglik=loglik, formula=formula, ode.opts=ode.opts))

        f.env <- new.env()
        ## set initial values
        assign("oldnll",NULL,f.env)
        assign("oldpar",NULL,f.env)
        assign("oldgrad",NULL,f.env)
        assign("data", data, f.env)

        objfun <- function(par, formula, model, loglik, observation, times, ode.opts) {
            if (identical(par,oldpar)) {
                if (debug) cat("returning old version of value\n")
                return(oldnll)
            }
            if (debug) cat("computing new version (nll)\n")

            v <- ode.sensitivity(par, formula, model, loglik, observation, times, ode.opts)
            oldnll <<- v[1]
            oldgrad <<- v[-1]
            oldpar <<- par

            return(oldnll)

        }
        gradfun <- function(par, formula, model, loglik, observation, times, ode.opts) {
            if (identical(par,oldpar)) {
                if (debug) cat("returning old version of grad\n")
                return(oldgrad)
            }
            if (debug) cat("computing new version (grad)\n")

            v <- ode.sensitivity(par, formula, model, loglik, observation, times, ode.opts)
            oldnll <<- v[1]
            oldgrad <<- v[-1]
            oldpar <<- par
            return(oldgrad)
        }

        environment(objfun) <- f.env

        environment(gradfun) <- f.env

        parnames <- c(newpar, loglik@par)
        attr(objfun, "parnames") <- parnames

        m <- mle2(objfun,
                  vecpar=TRUE,
                  start=start,
                  method="BFGS",
                  control=control,
                  gr=gradfun,
                  data=dataarg)

        .Object@mle2 <- m

        coef <- coef(m)

        if (!missing(links)) {
            coef <- transpar(coef, transform, inverse)
            thess <- numDeriv::jacobian(gradfun, coef, formula=orig.formula,
                model=orig.model,loglik=loglik,
                observation=data[,2],
                times=data[,1],
                ode.opts=ode.opts)
            vcov <- base::solve(thess)
            colnames(vcov) <- rownames(vcov) <- names(coef)
        } else {
            vcov <- vcov(m)
        }

        .Object@coef <- coef
        .Object@vcov <- vcov
        .Object@min <- m@min

        .Object
    }
)

##' Calculate sensitivity of the likelihood function with respect to the parameters
##'
##' @param parms ode/likelihood function parameters
##' @param formula formula specifing observation variable and mean
##' @export
ode.sensitivity <- function(parms, formula,
                        model, loglik,
                        observation, times=NULL,
                        ode.opts=list(method="rk4", hini=0.1)) {
    if (is.null(times)) times <- seq(length(observation))
    solution <- solve(model, times, parms, ode.opts=ode.opts)
    expr <- as.expression(formula[[3]])

    frame <- c(solution@solution, parms)

    mean <- eval(expr, frame)

    nstate <- length(model@state)

    dmds <- lapply(model@state, function(s) Deriv(expr, s))
    dmdp <- lapply(model@par, function(p) Deriv(expr, p))

    sens <- vector('list', nstate)
    for(i in 1:nstate) {
        sens[[i]] <- eval(dmds[[i]], frame) * solution@sensitivity[[i]]
    }

    sens_p <- sapply(dmdp, eval, frame)

    sens <- do.call("+", sens)

    if(class(sens_p) == "list")
        sens <- sens + do.call("cbind", sens_p)

    loglik.par <- as.list(parms[-c(1:length(model@par))])

    nll <- -sum(Eval(loglik, observation, mean, loglik.par))
    loglik.gr <- grad(loglik, observation, mean, loglik.par)
    sensitivity <- c(-colSums(loglik.gr[[1]] * sens))
    if(length(loglik.gr) > 1) sensitivity <- c(sensitivity, -sapply(loglik.gr[-1], sum))

    c(nll, sensitivity)
}

