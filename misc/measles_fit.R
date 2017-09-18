if(!require(fitode)) {
    devtools::install_github("parksw3/fitode"); library(fitode)
}

measles_df <- read.csv(url("https://ms.mcmaster.ca/~bolker/measdata/ewmeas.dat"), header=FALSE, sep=" ")
measles_df <- setNames(measles_df, c("year", "cases"))

SEI <- new("model.ode",
    name="SEI",
    model=list(
        S~mu*(N-S)-beta*I*S,
        E~beta*I*S-(mu+Sigma)*E,
        I~Sigma*E-(mu+gamma)*I
    ),
    initial=list(
        S~N*(1-e0-i0-r0),
        E~N*e0,
        I~N*i0
    ),
    par=c("N", "mu", "beta", "Sigma", "gamma", "e0", "i0", "r0")
)

measles <- Transform(
    SEI,
    transforms=list(
        beta~c0+c1*(1+cos(2*pi*t))
    ),
    par=c("N", "mu", "c0", "c1", "Sigma", "gamma", "e0", "i0", "r0")
)

start <- c(
    N=2e7,
    mu=0.02,
    c0=1.7e-5,
    c1=2e-5,
    Sigma=35.6,
    gamma=33,
    e0=1e-10,
    i0=1e-7,
    r0=0.94,
    sigma=9000
)

tdiff <- diff(measles_df$year)
tdiff <- c(tdiff, tail(tdiff,1))

ff <- fitode(cases~Sigma*E*tdiff*0.67,
    start=start,
    model=measles,
    loglik=select_model("gaussian"),
    data=measles_df,
    tcol="year",
    link=list(
        N="log",
        mu="log",
        c0="log",
        c1="log",
        Sigma="log",
        gamma="log",
        e0="logit",
        i0="logit",
        r0="logit",
        sigma="log"
    ),
    control=list(maxit=1e3),
    skip.hessian=TRUE,
    ode.opts=list(method="rk4"),
    debug=TRUE
)

start_nbinom <- coef(ff)[1:(length(coef(ff))-1)]
start_nbinom[["k"]] <- 10

fitode:::logLik.sensitivity(
    start_nbinom,
    cases~Sigma*E*1/52*0.67,
    measles,
    select_model("nbinom"),
    measles_df$cases,
    measles_df$year
)

## FIXME: return an error when link names does not match names of the parameters

ff_nbinom <- fitode(cases~Sigma*E*1/52*0.67,
    start=start_nbinom,
    model=measles,
    loglik=select_model("nbinom"),
    data=measles_df,
    tcol="year",
    link=list(
        N="log",
        mu="log",
        c0="log",
        c1="log",
        Sigma="log",
        gamma="log",
        e0="logit",
        i0="logit",
        r0="logit",
        k="log"
        ),
    control=list(maxit=1e3),
    skip.hessian=TRUE,
    ode.opts=list(method="rk4"),
    debug=TRUE
)

save("ff", "ff_nbinom", file="measles_fit.rda")

