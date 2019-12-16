library(fitode)

SI_model <- odemodel(
    name = "SI",
    model = list(
        S ~ - beta*S*I/N,
        I ~ beta*S*I/N - gamma*I
    ),
    observation = list(
        infected ~ ols(mean=I)
    ),
    initial = list(
        S ~ N * (1 - i0),
        I ~ N * i0
    ),
    par=c("beta", "gamma", "N", "i0"),
    link=c(i0="logit")
)

parms <- c(beta=1,gamma=0.5, N=1000, i0=1e-2, sd=10)

df <- simulate(SI_model, times=1:20, parms=parms, seed=101)

df$infected <- df$I + rnorm(length(df$I), sd=10)

ff <- fitode(
    model=SI_model,
    data=df,
    start=parms,
    fixed=c(N=1000)
)

ff2 <- update(ff,
    observation = list(
        log(infected) ~ dnorm(mean=log(I), sd=sd)
    ),
    par=c(SI_model@par, "sd")
)

ff
ff2

ff@vcov
ff2@vcov
