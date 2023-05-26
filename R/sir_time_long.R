library(fitode)
set.seed(101)

## hard-coded parameters to illustrate exogenous birth forcing; these probably can't be optimized
## functions involving parameters for optimization need to be inline in the formula
## could also convert t to an index for lookup in an existing vector
## (or interpolate, or whatever)
pc_birth_fun <- function(t) {
    return(0.2 + ((2050-t)/200)^2)
}

## model with exogenous per capita birth rate, sinusoidally/seasonally varying contact rate
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

## reasonable parameters (taken from example in epidemic/fitode vignettes for
## Sierra Leone Ebola example)
SIR_start <- c(beta0=70, beta1 = 0.1, gamma=60, N=40000, i0=0.0004, phi=6)

## biweekly 'data' for 200 years
system.time(ss_SIR <- simulate(SIR_model,
                               parms=SIR_start, times=seq(2014, 2050, by = 1/26)))
## ~1 seconds to simulate

## solving at a weekly time scale to improve numerical stability
system.time(ss_SIR2 <- simulate(SIR_model,
                               parms=SIR_start, times=seq(2014, 2050, by = 1/26),
                               solver.opts=list(method="rk4", hini=1/52)))
## ~ 2 seconds

## daily
system.time(ss_SIR3 <- simulate(SIR_model,
                                parms=SIR_start, times=seq(2014, 2050, by = 1/26),
                                solver.opts=list(method="rk4", hini=1/365)))
## ~ 11 seconds

# they all predict slightly different dynamics... not surprising!
plot(R~times, data=ss_SIR, type="l", log="y")
lines(R~times, data=ss_SIR2, col=2)
lines(R~times, data=ss_SIR3, col=3)

## for improved stability?
SIR_log_model <- odemodel(
    name="SIR (nbinom)",
    model=list(
        S ~ pc_birth_fun(t)*N - beta0*(1+beta1*cos(2*pi*t)) * S * exp(logI)/N,
        logI ~ beta0*(1+beta1*cos(2*pi*t)) * S/N - gamma,
        R ~ gamma * exp(logI)
    ),
    observation=list(
        confirmed ~ dnbinom(mu=R, size=phi)
    ),
    initial=list(
        S ~ N * (1 - i0),
        logI ~ log(N * i0),
        R ~ 0
    ),
    diffnames="R",
    par=c("beta0", "beta1", "gamma", "N", "i0", "phi"),
    link=c(i0="logit", beta0 = "log", beta1 = "logit", gamma = "log")
)

## biweekly 'data'
system.time(ss_log_SIR <- simulate(SIR_log_model,
                               parms=SIR_start, times=seq(2014, 2050, by = 1/26)))

plot(R~times, data=ss_SIR2, type="l", log="y") ## weekly simulation without logging
lines(R~times, data=ss_log_SIR, type="l", log="y", col=2) ## biweekly simulation with logging

## biweekly 'data' for 200 years
system.time(ss_log_SIR_long <- simulate(SIR_log_model,
                                   parms=SIR_start, times=seq(2014, 2214, by = 1/26)))

plot(R~times, data=ss_log_SIR_long, type="l", log="y")

## this somehow works although it gives a really bad fit.
system.time(SIRfit <- fitode(SIR_log_model, data=ss_log_SIR_long,
                             start=SIR_start,
                             control = list(maxit = 1e5),
                             tcol="times"))
## 600 seconds

plot(SIRfit)

sum(dnbinom(ss_log_SIR_long$confirmed, mu=ss_log_SIR_long$R, size=coef(SIRfit)[["phi"]], log=TRUE), na.rm=TRUE)
