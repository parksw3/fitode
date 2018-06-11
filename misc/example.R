library(fitode)

SI_model <- new("model.ode",
    name = "SI",
    model = list(
        S ~ - beta*S*I/N,
        I ~ beta*S*I/N - gamma*I
    ),
    observation = list(
        susceptible ~ dnorm(mean=S, sd=sigma1),
        infected ~ dnorm(mean=I, sd=sigma2)
    ),
    initial = list(
        S ~ N * (1 - i0),
        I ~ N * i0
    ),
    par=c("beta", "gamma", "N", "i0", "sigma1", "sigma2")
)

tvec <- 1:40

f <- ode.solve(SI_model, tvec, c(beta=1,gamma=0.5, N=100, i0=1e-3, sigma1=0.1, sigam2=0.1))

set.seed(101)
df <- data.frame(
    day=tvec,
    susceptible=rnorm(length(tvec), mean=f@solution$S, sd=0.1),
    infected=rnorm(length(tvec), mean=f@solution$I, sd=0.1)
)

start <- c(beta=1, gamma=0.5, N=100, i0=1e-3, sigma1=0.1, sigma2=0.1)

system.time(ff <- fitode(
    model=SI_model,
    data=df,
    start=start,
    tcol="day",
    link = list(
        beta="log",
        gamma="log",
        N="log",
        i0="logit",
        sigma1="log",
        sigma2="log"
    )
))

system.time(pp <- profile(ff@mle2, continuation="naive", trace=TRUE))

SI_model_R0 <- Transform(
    SI_model,
    list(beta~R0*gamma),
    par=c("R0", "gamma", "N", "i0")
)

SI_model_R0_u <- Transform(
    SI_model_R0,
    list(R0~1+RR),
    par=c("RR", "gamma", "N", "i0")
)

ff2 <- fitode(Deaths~gamma*I,
    start=c(RR=1, gamma=1, N=1e5, i0=1e-4, k=5),
    model=SI_model_R0_u,
    loglik=select_model("nbinom"),
    data=harbin,
    tcol="week",
    links = list(
        RR="log",
        gamma="log",
        N="log",
        i0="logit"
    )
)

bombay <- fitsir::bombay
bombay2 <- rbind(bombay, data.frame(week=32, mort=NA))

SI_model_c <- new("model.ode",
    name = "SI",
    model = list(
        S ~ - beta*S*I/N,
        I ~ beta*S*I/N - gamma*I,
        cDeath ~ gamma*I
    ),
    initial = list(
        S ~ N * (1 - i0),
        I ~ N * i0,
        cDeath ~ 0
    ),
    par=c("beta", "gamma", "N", "i0")
)

start <- c(beta=2, gamma=1, N=3000, i0=1e-5, ll.sigma=1)

system.time(ff4 <- fitode(mort|week ~ .diff(cDeath),
    start=start,
    model=SI_model_c,
    loglik=select_model("gaussian"),
    data=bombay2,
    link = list(
        beta="log",
        gamma="log",
        N="log",
        i0="logit"
    )
))
