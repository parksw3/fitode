library(fitode)

harbin <- fitsir::harbin

SI_model <- new("model.ode",
    name = "SI",
    model = list(
        S ~ - beta*S*I/N,
        I ~ beta*S*I/N - gamma*I
    ),
    initial = list(
        S ~ N * (1 - i0),
        I ~ N * i0
    ),
    par=c("beta", "gamma", "N", "i0")
)

start <- c(beta=2, gamma=1, N=1e5, i0=1e-4, sigma=5)

system.time(ff <- fitode(Deaths|week~gamma*I,
    start=start,
    model=SI_model,
    loglik=select_model("gaussian"),
    data=harbin,
    link = list(
        beta="log",
        gamma="log",
        N="log",
        i0="logit"
    )
))

system.time(ff2 <- fitsir::fitsir(
    data=harbin,
    start=start,
    type="death",
    icol="Deaths",
    tcol="week",
    optimizer="optim",
    method="BFGS"
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

harbin2 <- rbind(harbin, data.frame(week=19, Deaths=NA))

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

system.time(ff4 <- fitode(Deaths|week ~ .diff(cDeath),
    start=start,
    model=SI_model_c,
    loglik=select_model("gaussian"),
    data=harbin2,
    link = list(
        beta="log",
        gamma="log",
        N="log",
        i0="logit"
    )
))

