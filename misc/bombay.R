## exploring Bombay-fitting problems

## 1. Recapitulate fitting procedures from the paper
library(fitode)
library(ggplot2); theme_set(theme_bw())

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

bombay2 <- rbind(
    c(times=bombay$week[1] -
          diff(bombay$week)[1], mort=NA),
    bombay
)

## SIR_start <- c(beta=beta.nls, gamma=gamma.nls, I0=I0.KM, S0=S0.nls, k=50)
SIR_start <- c(beta = 7.89e-05, gamma = 4.47, I0 = 1, S0 = 61500, k = 50)
SIR_fit <- fitode(
    model = SIR_model,
    data = bombay2,
    start = SIR_start,
    fixed = list(gamma=SIR_start[["gamma"]]),
    tcol = "week",
    control = list(parscale = SIR_start[setdiff(names(SIR_start), "gamma")])
)

coef.bombay.fitode.nb <- coef(SIR_fit)
ci.bombay.fitode.nb <- confint(SIR_fit)[,-1]

gamma.nls <- SIR_start[["gamma"]]
N.KM <- 10^6
ci.derived.bombay.fitode.nb <- confint(SIR_fit,
        parm=list(
            Reff~beta*S0/gamma.nls,
            R0~beta*N.KM/gamma.nls,
            Tg~1/gamma.nls
        ))

cov2cor(vcov(SIR_fit))

## implausible CIs
ci.derived.bombay.fitode.nb

tmax <- 33
tvals <- 0:tmax
SIR_confband <- predict(SIR_fit, times=tvals, level=0.95)$mort

setup_plot <- function(xlab="", ylab="deaths",
                       at=100*(0:10), ...) {
    plot(NA, NA , bty="L", type="n",
         xaxs="i", yaxs="i", las=1,
         yaxt="n",
         xlab=xlab, ylab=ylab,
         ...
         )
    axis(side=2, at=at, las=1)
}
draw_confband <- function(confband, col = col.fitodenbCI) {
    with(confband,{
        polygon(
            x = c(times, rev(times)),
            y = c(`2.5 %`, rev(`97.5 %`)),
            col = col,
            border = NA,
            xpd = NA
        )
    })
}

col.fitodenbCI <- adjustcolor("black", alpha = 0.2)

tvals <- seq(0, tmax, length.out=1000)
setup_plot(ylim = c(0, 1500), xlim=c(0,tmax),
           xlab="weeks", ylab="plague deaths")
with(SIR_confband, lines(times, estimate))
with(bombay2, points(week, mort))

## yikes!
draw_confband(SIR_confband, col=col.fitodenbCI)

## Steve Ellner gives this as an exercise here:
##  https://ms.mcmaster.ca/~bolker/eeid/2010/Ecology/EllnerMaxLikEEID2010.pdf

## desperation: fix gamma, I0, *and* S0
SIR_fit2 <- fitode(
    model = SIR_model,
    data = bombay2,
    start = SIR_start,
    fixed = list(gamma=SIR_start[["gamma"]], I0 = 1, S0 = 6.2e4),
    tcol = "week",
    control = list(parscale = SIR_start[setdiff(names(SIR_start), c("gamma", "I0", "S0"))])
)

confint(SIR_fit,
        parm=list(
            Reff~beta*S0/gamma.nls,
            R0~beta*N.KM/gamma.nls,
            Tg~1/gamma.nls
        ))
coef(SIR_fit2)
cov2cor(vcov(SIR_fit2))


## Trying brute-force/from-scratch: make sure problems aren't somehow fitode bugs ...
sir.grad <- function (t, x, params) {
    with(as.list(c(x, params)), {
         dS.dt <- -beta*S*I
         dI.dt <- beta*S*I-gamma*I
         dR.dt <- gamma*I
         dxdt <- c(dS.dt,dI.dt,dR.dt) 
         list(dxdt)
    })
}

sir.pred <- function(beta, gamma, S0, I0) {
    times <- bombay2$week
    ode.params <- c(beta=beta, gamma=gamma)
    xstart <- c(S=S0,I=I0,R=0)
    out <- deSolve::ode(func=sir.grad,
               y=xstart,
               times=times,
               parms=ode.params
               ) |> as.data.frame()
    out
}

sir.nll <- function(beta, gamma, S0, I0, k) {
    out <- sir.pred(beta, gamma, S0, I0)
    ## not sure about starting value etc
    -sum(dnbinom(x = bombay2$mort[-1], mu = diff(out$R)[-1], size = k, log = TRUE))
}

fixvars <- c("gamma", "I0")
do.call(sir.nll, as.list(SIR_start))
m1 <- mle2(sir.nll, start = as.list(SIR_start), method = "Nelder-Mead",
           control = list(parscale = abs(SIR_start)[setdiff(names(SIR_start), fixvars)], maxit = 1e5),
           fixed = as.list(SIR_start)[fixvars]
           )
cov2cor(vcov(m1))

pred <- with(as.list(coef(m1)), sir.pred(beta, gamma, S0, I0))
plot(diff(pred$R), type ="l")
with(bombay2, points(week, mort))
with(predict(SIR_fit2)$mort, lines(times, estimate, col = 2))
with(predict(SIR_fit)$mort, lines(times, estimate, col = 4))

rng <- 0.2
liksurf <- expand.grid(beta = seq(1-rng, 1+rng, length = 51)*coef(m1)[["beta"]],
                       S0   = seq(1-rng, 1+rng, length = 51)*coef(m1)[["S0"]])

## vary only beta and S0
liksurf$lik <- apply(liksurf, 1,
                     function(p) {
                         sir.nll(beta = p[[1]], gamma = coef(m1)[["gamma"]],
                                 S0 = p[[2]], I0 = coef(m1)[["I0"]],
                                 k = coef(m1)[["k"]])
                     })

liksurf$slik <- liksurf$lik-min(liksurf$lik, na.rm= TRUE)

ggplot(liksurf, aes(beta, S0)) +
    geom_raster(aes(fill = slik)) +
    scale_fill_continuous(trans = "log10") +
    geom_contour(aes(z = slik), breaks  = 2, colour = "red")

pp <- profile(m1)
confint(pp)

ppb <- as.data.frame(pp@profile$beta)


## S0 is basically == N, so strongly confounded with beta

fixvars <- c("S0", "I0")
m2 <- mle2(sir.nll, start = as.list(SIR_start), method = "Nelder-Mead",
           control = list(parscale = abs(SIR_start)[setdiff(names(SIR_start), fixvars)], maxit = 1e5),
           fixed = as.list(SIR_start)[fixvars]
           )

rng <- 0.2
liksurf <- expand.grid(beta = seq(1-rng, 1+rng, length = 51)*coef(m2)[["beta"]],
                       gamma   = seq(1-rng, 1+rng, length = 51)*coef(m2)[["gamma"]])

## vary only beta and gamma
liksurf$lik <- apply(liksurf, 1,
                     function(p) {
                         sir.nll(beta = p[[1]],
                                 gamma = p[[2]],
                                 S0 = coef(m2)[["S0"]],
                                 I0 = coef(m2)[["I0"]],
                                 k = coef(m1)[["k"]])
                     })

liksurf$slik <- liksurf$lik-min(liksurf$lik, na.rm= TRUE)

ggplot(liksurf, aes(beta, gamma)) +
    geom_raster(aes(fill = slik)) +
    scale_fill_continuous(trans = "log10") +
    geom_contour(aes(z = slik), breaks  = 2, colour = "red")
