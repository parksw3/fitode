if(!require(fitode)) {
    devtools::install_github("parksw3/fitode"); library(fitode)
}

measles_df <- read.csv(url("https://ms.mcmaster.ca/~bolker/measdata/ewmeas.dat"), header=FALSE, sep=" ")
measles_df <- setNames(measles_df, c("year", "cases"))

SEI_obs <- new("model.ode",
               name="SEI observe",
               model=list(
                   S~mu*(N-S)-beta*(I.obs+I.unobs)*S,
                   E~beta*(I.obs+I.unobs)*S-(mu+sigma)*E,
                   I.unobs~sigma*E-(mu+gamma+lambda)*I.unobs,
                   I.obs~lambda*I.unobs-(mu+gamma)*I.obs
               ),
               initial=list(
                   S~N*(1-e0-i0.unobs-i0.obs-r0),
                   E~N*e0,
                   I.unobs~N*i0.unobs,
                   I.obs~N*i0.obs
               ),
               par=c("N", "mu", "beta", "sigma", "gamma", "lambda", "e0", "i0.unobs", "i0.obs", "r0")
)

measles_obs <- Transform(
    SEI_obs,
    transforms=list(
        beta~c0+c1*(1+cos(2*pi*(t-t0)))
    ),
    par=c("N", "mu", "c0", "c1", "sigma", "gamma", "lambda", "t0", "e0", "i0.unobs", "i0.obs", "r0")
)

start <- c(
    N=5e7,
    mu=0.02,
    c0=1.7e-5,
    c1=3e-5,
    sigma=35.6,
    gamma=33,
    lambda=10,
    t0=0.01,
    e0=1e-10,
    i0.unobs=1e-7,
    i0.obs=1e-8,
    r0=0.94,
    sigma=10
)

ff <- fitode(cases~lambda*I.obs*1/52,
             start=start,
             model=measles_obs,
             loglik=select_model("gaussian"),
             data=measles_df,
             tcol="year",
             link=list(
                 N="log",
                 mu="log",
                 c0="log",
                 c1="log",
                 sigma="log",
                 gamma="log",
                 lambda="log",
                 t0="logit",
                 e0="logit",
                 i0.unobs="logit",
                 i0.obs="logit",
                 r0="logit",
                 sigma="log"
             ),
             control=list(maxit=1e3),
             skip.hessian=TRUE
)

save("ff", file="measles_fit.rda")

