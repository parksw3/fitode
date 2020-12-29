##' @importFrom MASS mvrnorm
##' @importFrom mvtnorm dmvnorm
##' @importFrom stats approx
##' @keywords internal
wmvrnorm <- function(object,
                     nsim=1000,
                     sample.vcov,
                     seed) {
    if(!missing(seed)) set.seed(seed)

    model <- Transform(object@model, keep_sensitivity=FALSE)

    vv <- vcov(object, "links")
    vv[lower.tri(vv)] <- t(vv)[lower.tri(vv)]

    simtraj <- vector('list', length(model@expr))

    for (j in 1:length(model@expr)) simtraj[[j]] <- matrix(NA, nrow=length(unique(object@data$times)), ncol=nsim)

    traj.logLik <- rep(NA, nsim)

    simpars <- MASS::mvrnorm(nsim,mu=coef(object, "links"),
                             Sigma=vv)

    sample.logLik <- mvtnorm::dmvnorm(simpars, coef(object, "links"), vv, log=TRUE)

    simpars_orig <- t(apply(simpars, 1, apply_link, object@mle2@data$linklist, "linkinv"))

    for (i in 1:nsim) {
        ss.tmp <- logLik.sensitivity(simpars_orig[i,], model, object@data,
                                     solver.opts=object@mle2@data$solver.opts,
                                     solver=object@mle2@data$solver,
                                     return.traj=TRUE)
        for (k in 1:length(model@expr)) {
            simtraj[[k]][,i] <- ss.tmp$traj[[k]]
        }

        if (length(object@prior) > 0) {
            lprior <- eval(object@mle2@data$priorlist$prior.density, as.list(simpars[i,]))
        } else {
            lprior <- 0
        }

        traj.logLik[i] <- -ss.tmp$nll[1] + lprior
    }
    log.ww <- traj.logLik-sample.logLik

    ## FIXME: prevent taking exponential causing underflow/overflow (hopefully??)
    ww <- exp(log.ww-median(log.ww))
    ww <- ww/sum(ww)

    list(
        simtraj=simtraj,
        simpars_orig=simpars_orig,
        weight=ww
    )
}

## wquant from King et al.
wquant <- function (x, weights, probs = c(0.025, 0.975)) {
    which <- !is.na(weights)
    x <- x[which]
    weights <- weights[which]

    if (all(is.na(x)) || length(x) == 0) return(rep(NA, length(probs)))

    idx <- order(x)
    x <- x[idx]
    weights <- weights[idx]
    w <- cumsum(weights)/sum(weights)
    rval <- approx(w,x,probs,rule=1)
    rval$y
}
