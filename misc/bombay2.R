## general utils
xfun <- function(x, n) x[setdiff(names(x), n)]

## gradient function
options(sir.grad_param = "R0m1_r")
sir.grad <- function (t, x, params) {
    with(as.list(c(x, params)), {
        parameterization <- getOption("sir.grad_param", "R0m1_r")
        if (parameterization == "R0m1_r") {
            R0m1 <- exp(logR0m1)
            r <- exp(logr)
            gamma <- r/R0m1
            beta <- (r + gamma)/N
        } else if (parameterization == "R0m1_gamma") {
            R0m1 <- exp(logR0m1)
            gamma <- exp(loggamma)
            beta <- (R0m1+1)*gamma/N
        } else if (parameterization == "beta_gamma") {
            gamma <- exp(loggamma)
            beta <- exp(logbeta)
        }
        dS.dt <- -beta*S*I
        dI.dt <- beta*S*I-gamma*I
        dR.dt <- gamma*I
        dxdt <- c(dS.dt,dI.dt,dR.dt) 
        list(dxdt)
    })
}
sir.pred <- function(p) {
    times <- bombay2$week
    I0 <- exp(p[["logI0"]])
    S0 <- exp(p[["logS0"]])
    ode.params <- c(xfun(p, c("logI0", "logS0", "logk")), N = S0 + I0)
    xstart <- c(S=S0, I=I0, R=0)
    out <- deSolve::ode(func=sir.grad,
                        y=xstart,
                        times=times,
                        parms=ode.params
                        ) |> as.data.frame()
    diff(out$R)
}

library(deSolve)
data("bombay", package = "fitode")
## incidence-based version
bombay2 <- rbind(
    c(times=bombay$week[1] -
          diff(bombay$week)[1], mort=NA),
    bombay
)

## cheat a bit: use starting values from previous fit
##  don't know how badly we'll do if we start somewhere else
## (use DEoptim or ... ???)

## R0m1 parameterization
start_R0m1_gamma <- c(logR0m1 = -2.4342887050396,
            loggamma = 1.56234630490025, 
            logS0 = 10.939429567809, 
            logI0 = -0.19938133036307,
            logk = 3.91371837024452)

## no parameters constrained/fixed
sir_nll_all <- function(p, param = NULL) {
    if (!is.null(param)) {
        op <- options(sir.grad_param = param)
        on.exit(options(op))
    }
    sirpred <- sir.pred(xfun(p, c("logk")))
    -sum(dnbinom(x = bombay2$mort[-1],
                 mu = sirpred,
                 size = exp(p["logk"]), log = TRUE))
}

options(sir.grad_param = "R0m1_gamma")
sir_nll_all(start_R0m1_gamma)
opt_R0m1_gamma <- optim(start_R0m1_gamma,
              fn = sir_nll_all,
              control = list(maxit = 1e6),
              hessian = TRUE)

coefraw_R0m1_gamma  <- data.frame(est = opt_R0m1_gamma$par,
                   se = sqrt(diag(solve(opt_R0m1_gamma$hessian))))
print(coefraw_R0m1_gamma)

coefresp_R0m1_gamma <- transform(coefraw_R0m1_gamma,
                                       se = se*exp(est),
                                       est = exp(est))
print(coefresp_R0m1_gamma)

## R0m1_gamma <-> beta_gamma
trfun1 <- function(p, out = "beta") {
    logN <- log(exp(p[["logS0"]]) + exp(p[["logI0"]]))
    if (out == "beta") {
        beta_est <- with(as.list(p), {
            log1p(exp(logR0m1))+loggamma-logN
        })
        return(c(logbeta = beta_est, xfun(p, "logR0m1")))
    }
    if (out == "R0m1") {
        R0m1_est <- with(as.list(p), {
            log(exp(logbeta+logN-loggamma)-1)
        })
        return(c(logR0m1 = R0m1_est, xfun(p, "logbeta")))
    }
    stop("unknown 'out' value")
}

## exp + name mutation
trfun2 <- function(x) {
    setNames(exp(x), gsub("^log", "", names(x)))
}
## round-trip test ...
stopifnot(all.equal(
    start_R0m1_gamma,
    trfun1(trfun1(start_R0m1_gamma, "beta"), "R0m1")))

start_beta_gamma <- trfun1(opt_R0m1_gamma$par)

options(sir.grad_param = "beta_gamma")

## we get the same NLL from the translated 
stopifnot(all.equal(sir_nll_all(start_beta_gamma),
                    sir_nll_all(opt_R0m1_gamma$par, param = "R0m1_gamma")))

opt_beta_gamma <- optim(start_beta_gamma,
              fn = sir_nll_all,
              control = list(maxit = 1e6),
              hessian = FALSE)

## ?? degenerate Nelder-Mead simplex ???
## almost but not quite identical ... beta achieves slightly *lower* NLL
cbind(R0m1 = opt_R0m1_gamma$value,
      beta = opt_beta_gamma$value)

## lower nll but *nearly* identical params
cbind(start_beta_gamma, opt_beta_gamma$par)

## back-transform results to R0m1_gamma param
obg_R0m1 <- trfun1(opt_beta_gamma$par, out = "R0m1")
stopifnot(all.equal(opt_beta_gamma$value,
                    sir_nll_all(obg_R0m1, param = "R0m1_gamma")))

## scales are not crazy, so I don't think parscale will help ...
## (range(abs(x)) approx. 4-17)
try(optimHess(opt_beta_gamma$par, sir_nll_all))

## we can compute it with tiny ndeps but it's silly
## (similar results with smaller ndeps
try(optimHess(opt_beta_gamma$par, sir_nll_all,
              control = list(ndeps = rep(1e-4, 5))))
## doesn't throw an error, but bogus results ...
numDeriv::hessian(sir_nll_all, opt_beta_gamma$par)
## could fart around with method.args() but ugh ...

## but we **can** compute the Hessian in the R0m1_gamma space ...
optimHess(obg_R0m1, sir_nll_all, param = "R0m1_gamma")
## delta method to get back to beta/gamma scale...???

## perturb by hand?

## first check equal preds
options(sir.grad_param = "beta_gamma")
pred1 <- sir.pred(opt_beta_gamma$par)
options(sir.grad_param = "R0m1_gamma")
pred2 <- sir.pred(obg_R0m1)
stopifnot(all.equal(pred1, pred2)) ## of course ...

options(sir.grad_param = "beta_gamma")
## insane!!
pert1 <- opt_beta_gamma$par + c(0.01, rep(0,4))
trfun2(trfun1(pert1, out = "R0m1"))
trfun2(trfun1(opt_beta_gamma$par, out = "R0m1"))
## fairly radical difference implied in (R0-1) by a 0.01 change in logbeta ...

pred1B <- sir.pred(pert1)
options(sir.grad_param = "R0m1_gamma")
pred2B <- sir.pred(obg_R0m1 + + c(0.01, rep(0,4)))
all.equal(pred2, pred2B)  ## 2.7% difference

pertseq <- 10^seq(-5, -2, by = 0.5)
pertmat <- matrix(NA_real_, ncol = length(pertseq), nrow = length(pred1))
options(sir.grad_param = "beta_gamma")
for (i in seq_along(pertseq)) {
    pertmat[,i] <- sir.pred(opt_beta_gamma$par + c(pertseq[i], rep(0,4)))
}
par(las = 1, bty = "l")
matplot(pertmat, type = "b", pch=1, log = "y", ylim = c(0.1, 1e6),
        col = 2:8)
lines(pred1)
legend("topright",
       legend = c("orig",sprintf("10^(%1.1f)", log10(pertseq))),
       lty = 1,
       col = 1:8)


## plot some slices ... 4-dim (logS0, logI0, (logbeta/logR0m1), loggamma ...)
## re-inventing the slice() wheel from mle2 ...
calc_slice <- function(par, fun, rng = 0.1, nvec = NULL,
                       exclude = NULL, pb = TRUE, ...) {
    if (!is.null(exclude)) {
        xpars <- par[exclude]
        par <- xfun(par, names(xpars))
    }
    if (length(rng)<length(par)) rng <- rep(rng, length.out = length(par))
    if (is.null(nvec)) nvec <- c(rep(51, 2), rep(4, length(par)-2))
    parvecs <- Map(function(p, r, n) {
        seq(p*(1-r), p*(1+r), length.out = n)
    }, par, rng, nvec)
    parmat <- do.call(expand.grid, parvecs)
    ## progress bar?? parallelize?
    np <- nrow(parmat)
    if (pb) pbb <- txtProgressBar(max = np, style = 3)
    nll <- rep(NA_real_, np)
    for (i in seq_along(nll)) {
        if (pb) setTxtProgressBar(pbb, i)
        nll[i] <- fun(c(unlist(parmat[i,]), xpars), ...)
    }
    data.frame(parmat, nll)
}

## slow-ish (8-10 mins?)
fn <- "slices.rda"
if (file.exists(fn)) { load (fn) } else {
   system.time(cc1 <- calc_slice(obg_R0m1, exclude = "logk", sir_nll_all, param = "R0m1_gamma"))
   system.time(cc2 <- calc_slice(opt_beta_gamma$par, exclude = "logk", sir_nll_all, param = "beta_gamma"))
   ## try tiny range, coarser grid (for speed)
   system.time(cc3 <- calc_slice(opt_beta_gamma$par,
                                 rng = 0.0001,
                                 n = c(31, 31, 3, 3),
                                 exclude = "logk", sir_nll_all, param = "beta_gamma"))

   save(list = ls(pattern = "^cc[0-9]+"), file = fn)
                                 }


library(ggplot2); theme_set(theme_bw() + theme(panel.spacing=grid::unit(0, "lines")))
gg1 <- ggplot(cc1, aes(x=logR0m1, y = loggamma, fill = nll, z = nll)) +
    geom_raster() +
    facet_grid(logS0 ~ logI0) +
    scale_fill_viridis_c(trans = "log10") +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    geom_contour(breaks = min(cc1$nll, na.rm = TRUE) + 10, colour = "red")


print(gg1)

print(gg1 %+% cc2 + aes(x = logbeta))

## apparently we can only perturb the parameters a *tiny* bit and still get sensible answers ...
gg3 <- gg1 %+% cc3 + aes(x = logbeta)
## lots of warnings about dropping 'fill', I think because we have panels
##  that are all NA ...

suppressWarnings(print(gg3))

