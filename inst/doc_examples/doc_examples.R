library(fitode)
SIR_mort_model <- odemodel(
    name="SIR mortality rate",
    model=list(
        S ~ - beta * S * I/N,
        I ~ beta * S * I/N - gamma * I,
        R ~ gamma * I
    ),
    observation=list(
        mort ~ ols(mean = gamma * I)
    ),
    initial=list(
        S ~ N * (1 - i0),
        I ~ N * i0,
        R ~ 0
    ),
    par=c("beta", "gamma", "N", "i0"),
    link=c(i0="logit")
)
SIR_mort_start <- c(beta=0.8, gamma=0.4, N=6000, i0=0.001)
SIR_mort_fit <- fitode(
    SIR_mort_model,
    data=bombay,
    start=SIR_mort_start,
    tcol="week"
)
save(SIR_mort_fit, SIR_mort_start, SIR_mort_model,
     file = "SIR_mort.rda", version = 2)
