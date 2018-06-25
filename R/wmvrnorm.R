wmvrnorm <- function(object,
                     nsim=1000,
                     seed) {
    if(!missing(seed)) set.seed(seed)

    model <- Transform(object@model, keep_sensitivity=FALSE)

    vv <- vcov(object, "fitted")
    vv[lower.tri(vv)] <- t(vv)[lower.tri(vv)]

    simtraj <- vector('list', length(model@expr))

    for (j in 1:length(model@expr)) simtraj[[j]] <- matrix(NA, nrow=length(unique(object@data$times)), ncol=nsim)

    traj.logLik <- rep(NA, nsim)

    simpars <- MASS::mvrnorm(nsim,mu=coef(object, "fitted"),
                             Sigma=vv)

    sample.logLik <- mvtnorm::dmvnorm(simpars, coef(object, "fitted"), vv, log=TRUE)

    simpars_orig <- t(apply(simpars, 1, apply_link, object@mle2@data$linklist, "linkinv"))

    for (i in 1:nsim) {
        ss.tmp <- fitode:::logLik.sensitivity(simpars_orig[i,], model, object@data,
                                              solver.opts=object@mle2@data$solver.opts,
                                              solver=object@mle2@data$solver,
                                              fixed=object@fixed,
                                              return.traj=TRUE)
        for (k in 1:length(model@expr)) {
            simtraj[[k]][,i] <- ss.tmp$traj[[k]]
        }

        traj.logLik[i] <- -ss.tmp$nll[1]
    }
    log.ww <- traj.logLik-sample.logLik

    ## FIXME: prevent taking exponential causing underflow/overflow (hopefully??)
    ww <- exp(log.ww-median(log.ww))

    list(
        simtraj=simtraj,
        simpars_orig=simpars_orig,
        weight=ww
    )
}

## wquant from King et al.
wquant <- function (x, weights, probs = c(0.025, 0.975)) {
    idx <- order(x)
    x <- x[idx]
    weights <- weights[idx]
    w <- cumsum(weights)/sum(weights)
    rval <- approx(w,x,probs,rule=1)
    rval$y
}
