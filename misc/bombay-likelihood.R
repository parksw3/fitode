library(fitode)

bombay2 <- rbind(
    c(times=bombay$week[1] -
          diff(bombay$week)[1], mort=NA),
    bombay
)

SIR_model <- odemodel(
    name="SIR model",
    model=list(
        S ~ - beta * S * I,
        I ~ beta * S * I - gamma * I,
        R ~ gamma * I
    ),
    observation = list(
        mort ~ dnbinom(mu = R, size=k)
    ),
    diffnames="R",
    initial=list(
        S ~ S0,
        I ~ I0,
        R ~ 0
    ),
    par=c("beta", "gamma", "S0", "I0", "k")
)

SIR_fit2 <- fitode(
    SIR_model,
    data = bombay2,
    start = c(beta = 9.20444154491234e-05, S0 = 53113.5957792255, I0 = 0.870570340799881,
              k = 49.5911057134356),
    fixed = c(gamma = 1), # Tg = 1/gamma = 1 week
    tcol = "week"
)

nfits <- 41
gammavec <- seq(1, 14, length.out=nfits)
loglikvec <- rep(NA, nfits)
fitlist <- vector('list', length(gammavec))

for (i in 1:nfits) {
    if (i==1) {
        start_ll <- coef(SIR_fit2)
    } else {
        start_ll <- coef(fitlist[[which.max(loglikvec[1:(i-1)])]])
    }

    ## hini=1/4 to solve the ode at a finer scale.. so that things don't go
    ## crazy when we're starting with very high gamma
    fitlist[[i]] <- try(fitode(
        SIR_model,
        data = bombay2,
        start = start_ll,
        fixed = list(gamma=gammavec[i]),
        tcol = "week"
    ))

    if (!inherits(fitlist[[i]], "try-error")) {
        loglikvec[i] <- logLik(fitlist[[i]])
    }
}
