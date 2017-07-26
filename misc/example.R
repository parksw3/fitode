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
start <- c(log.beta=1.5, log.gamma=1, log.N=10, logit.i=-10, ll.k=2)

ff <- fitode(Deaths~I,
    start=start,
    model=SI_model_trans, loglik=select_model("nbinom"),
    data=harbin,
    tcol="week"
)

ff2 <- fitsir::fitsir(harbin, start, method="BFGS", dist="nbinom", tcol="week",icol="Deaths")

all.equal(coef(ff), coef(ff2), tolerance = 1e-3) ## returns FALSE if we lower the tolerance

matplot(harbin$week, predict(ff, level=0.95, method="wmvrnorm")[,-1], type="l", col=1)
matlines(harbin$week, predict(ff2,level=0.95)[,-1], col=2)
points(harbin)

## incidence fitting

ff3 <- fitode(Deaths~exp(log.beta)*S*I/exp(log.N),
    start=start,
    model=SI_model_trans, loglik=select_model("nbinom"),
    data=harbin,
    tcol="week"
)

ff4 <- fitsir::fitsir(harbin, start, method="BFGS", dist="nbinom", tcol="week",icol="Deaths", type="incidence")

## we need high tolerance because fitode fits based on instantaneous incidence
## whereas fitsir fits based on actual incidence (difference in number of susceptible in consecutive observations)
all.equal(coef(ff3), coef(ff4), tolerance = 3e-2)

## difference is much smaller on a constrained scale
all.equal(
    fitsir::trans.pars(coef(ff3)),
    fitsir::trans.pars(coef(ff4)),
    tolerance=3e-3
)

matplot(harbin$week, predict(ff3, level=0.95)[,-1], type="l", col=1)
matlines(harbin$week, predict(ff4,level=0.95)[,-1], col=2)
points(harbin)

