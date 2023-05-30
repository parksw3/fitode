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
    link=c(i0="logit", beta1 = "logit")
)

SIR_model_t <- Transform(SIR_model,
          transform=list(beta0~R0*gamma),
          par=c("R0", "beta1", "gamma", "N", "i0", "phi"),
          link=c(i0="logit", beta1 = "logit"))

SIR_model_t@grad
