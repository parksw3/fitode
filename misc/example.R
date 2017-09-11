library(fitode)

harbin <- structure(list(week = 2:18,
                         Deaths = c(2, 7, 2, 6, 12, 68, 91,
                                    126, 229, 250, 212, 161, 101, 108, 46, 40, 14)),
                    .Names = c("week", "Deaths"),
                    row.names = c(NA, -17L), class = "data.frame")

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

start <- c(beta=2, gamma=1, N=1e5, i0=1e-4, k=5)

ff <- fitode(Deaths~gamma*I,
    start=start,
    model=SI_model,
    loglik=select_model("nbinom"),
    data=harbin,
    tcol="week",
    links = list(
        beta="log",
        gamma="log",
        N="log",
        i0="logit"
    )
)

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



