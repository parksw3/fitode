library(fitode)

## beta/gamma = R0
## beta = (R0_1+1)*gamma

SIR_model <- odemodel(
    name="SIR model",
    model=list(
        S ~ - beta * S * I,
        I ~ beta * S * I - gamma * I,
        R ~ gamma * I
    ),
    observation = list(
        mort ~ dnbinom(mu = R, size = k)
    ),
    diffnames="R",
    initial=list(
        S ~ S0,
        I ~ I0,
        R ~ 0
    ),
    par=c("beta", "gamma", "S0", "I0", "k")
)

SIR_model_reparam <- Transform(SIR_model,
          transform=list(beta~(R0_1+1)*gamma/(S0+I0)),
          par=c("R0_1", "gamma", "S0", "I0", "k"),
          link=c(R0_1="log"))

bombay2 <- rbind(
    c(times=bombay$week[1] -
          diff(bombay$week)[1], mort=NA),
    bombay
)

SIR_start <- c(R0_1=3.02e-4, gamma=1.35e3, I0=3.11e-3, S0=1.46e7, k=58.4)

SIR_start2 <- c(R0_1=3.02e-5, gamma=1.35e3, I0=3.11e-3, S0=1.46e8, k=58.4)

SIR_start3 <- c(R0_1=8.953650e-05, gamma=4.560758e+03, I0=9.206802e-04, S0=4.942489e+07, k=60)

ss <- simulate(SIR_model_reparam, par=SIR_start3, times=bombay2$week)

plot(ss$R, type="l")

SIR_fit <- fitode(
    model = SIR_model_reparam,
    data = bombay2,
    start = SIR_start,
    tcol = "week"
)

SIR_fit2 <- fitode(
    model = SIR_model_reparam,
    data = bombay2,
    start = SIR_start2,
    tcol = "week"
)

SIR_fit3 <- fitode(
    model = SIR_model_reparam,
    data = bombay2,
    start = SIR_start3,
    tcol = "week"
)
