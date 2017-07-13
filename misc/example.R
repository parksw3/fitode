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
    par= c("beta", "gamma", "N", "i0")
)
SI_model_trans <- Transform(SI_model,
    transforms=list(
        beta~exp(log.beta),
        gamma~exp(log.gamma),
        N~exp(log.N),
        i0~(1+tanh(logit.i/2))/2
    ),
    par=c("log.beta", "log.gamma", "log.N","logit.i")
)
start <- c(log.beta=1.5, log.gamma=1, log.N=10, logit.i=-10)

ff <- fitode(Deaths~I,
    start=start,
    model=SI_model_trans, loglik=select_model("poisson"),
    data=harbin,
    tcol="week"
)

ff2 <- fitsir::fitsir(harbin, start, method="BFGS", dist="poisson", tcol="week",icol="Deaths")

all.equal(coef(ff), coef(ff2), tolerance = 1e-4) ## returns FALSE if we lower the tolerance

ss <- solve(SI_model_trans, harbin$week, coef(ff))

plot(harbin)
lines(ss@solution[,c(1,3)], type="l")
plot(ff2, add=TRUE, col.traj=2)
