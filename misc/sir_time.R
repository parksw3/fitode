library(fitode)
set.seed(101)

## hard-coded parameters to illustrate exogenous birth forcing; these probably can't be optimized
## functions involving parameters for optimization need to be inline in the formula
## could also convert t to an index for lookup in an existing vector
## (or interpolate, or whatever)
pc_birth_fun <- function(t) {
    return(0.2 + ((2050-t)/200)^2)
}

SIR_model <- odemodel(
    name="SIR (nbinom)",
    model=list(
        S ~ pc_birth_fun(t)*N - beta0*(1+beta1*cos(2*pi*t)) * S * I/N,
        I ~ beta0*(1+beta1*cos(2*pi*t)) * S * I/N - gamma * I,
        R ~ gamma * I
    ),
    observation=list(
        confirmed ~ dnbinom(mu=R, size=phi)
    ),
    initial=list(
        S ~ N * (1 - i0),
        I ~ N * i0,
        R ~ 0
    ),
    diffnames="R",
    par=c("beta0", "beta1", "gamma", "N", "i0", "phi"),
    link=c(i0="logit", beta0 = "log", beta1 = "logit", gamma = "log")
)

SIR_start <- c(beta0=70, beta1 = 0.1, gamma=60, N=40000, i0=0.0004, phi=6)

## devtools::load_all("~/R/pkgs/fitode")
## biweekly 'data' for 200 years
system.time(ss_SIR <- simulate(SIR_model,
                               parms=SIR_start, times=seq(2014, 2214, by = 1/26)))
## 3.8 seconds

oss_sub <- ss_SIR[, setdiff(colnames(ss_SIR), c("times", "sim"))]
par(las = 1, bty = "l")
matplot(ss_SIR[,1], ss_sub, log = "y",type = "l", ylim = c(1e-1,1e6),
        ylab = "")
legend("top", horiz= TRUE,
       col = 1:(ncol(ss_sub)),
       lty = 1,
       legend = colnames(ss_sub)
       )

## debug(fitode)
system.time(SIRfit <- fitode(SIR_model, data=ss_SIR,
                             start=SIR_start,
                             control = list(maxit = 1e5),
                             trace = 100,
                             tcol="times"))
## sensitivities are *huge*

if (FALSE) {
    ## debugging
    debug(objgrad_fun)
    system.time(with(dataarg, objfun(start, data, solver.opts, solver, linklist, priorlist)))
    system.time(with(dataarg, gradfun(start, data, solver.opts, solver, linklist, priorlist)))
}
## fitode does **not** work with vector parameters at present!
if (FALSE) {
    SIR_model_vec <- odemodel(
    name="SIR (nbinom)",
    model=list(
        S ~ pc_birth_fun(t)*N - ifelse(t<2030, beta[1], beta[2])*(1+0.1*cos(2*pi*t)) * S * I/N,
        I ~ beta*(1+0.1*cos(2*pi*t)) * S * I/N - gamma * I,
        R ~ gamma * I
    ),
    observation=list(
        confirmed ~ dnbinom(mu=R, size=phi)
    ),
    initial=list(
        S ~ N * (1 - i0),
        I ~ N * i0,
        R ~ 0
    ),
    diffnames="R",
    par=c("beta", "gamma", "N", "i0", "phi", "nu"),
    link=c(i0="logit")
)

SIR_start_vec <- list(beta=c(70,60), gamma=60, N=40000, i0=0.0004, phi=6)

system.time(ss_SIR_vec <- simulate(SIR_model_vec,
    parms=SIR_start_vec, times=seq(2014, 2214, by = 1/26)))
}
