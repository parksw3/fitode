## ----setup, echo=FALSE, include=FALSE-----------------------------------------
## set default base-graphics plotting options
knitr::knit_hooks$set(basefig=function(before, options, envir) {
                   if (before) {
                       oldpar <- par(bty="l",las=1)
                       on.exit(par(oldpar))
                   } else { }
})
knitr::opts_chunk$set(error = FALSE)
## B&W format: https://stackoverflow.com/questions/22666475/how-to-achieve-blackwhite-friendly-syntax-highlighting-in-pdfs-produced-with-kn
##knitr::opts_knit$set( out.format="latex" )
##knitr::knit_theme$set("print")
## override BMB's finicky settings
options(warnPartialMatchArgs = FALSE,
          warnPartialMatchAttr = FALSE,
          warnPartialMatchDollar = FALSE)



## ----functions to compute KM parameters, echo=FALSE---------------------------
omega_fun <- function(Reff, gamma, S0, I0)
    {(gamma/2)*sqrt((Reff-1)^2+(2*I0/S0)*Reff^2)}
phi_fun <- function(Reff, gamma, S0, I0)
    {atanh((Reff-1)/(2*omega_fun(Reff,gamma,S0,I0)/gamma))}
a_fun <- function(Reff, gamma, S0, I0)
    {(2 * S0/(gamma * Reff^2))*omega_fun(Reff,gamma,S0,I0)^2}


## ----nls-bombay---------------------------------------------------------------
sech <- function(x) {1/cosh(x)}
KM_approx <- function(t, a, omega, phi) {a * sech(omega*t - phi)^2}
KM.parameters <- c(a = 890, omega = 0.2, phi = 3.4)
nlsfit <- nls(mort ~ KM_approx(week, a, omega, phi),
              data = fitode::bombay,
              start = KM.parameters)
nls.parameters <- coef(nlsfit)
print(nls.parameters)


## ----nls-bombay-2-------------------------------------------------------------
a.guess <- 1000    # crude "by eye" estimate of peak value,
tpeak.guess <- 15  # peak time,
omega.guess <- 1   # and half the initial growth rate
phi.guess <- omega.guess * tpeak.guess
nlsfit <- nls(mort ~ KM_approx(week, a, omega, phi),
              data = fitode::bombay,
              start = c(a = a.guess, omega = omega.guess,
                        phi = phi.guess))
print(nls.parameters <- coef(nlsfit))


## ----eval=FALSE---------------------------------------------------------------
# lower = c(a = 0, omega = 0, phi = 0)


## ----KM-fitted-values, echo=FALSE---------------------------------------------
dig <- 3 # number of significant digits for printing
assignfun <- function(par,
                      text=".KM") {
    parnames <- names(par)

    for (p in parnames) {
        assign(paste0(p, text), signif(par[[p]], dig),
               envir=globalenv())
    }
}

assignfun(KM.parameters)
tpeak.KM <- with(as.list(KM.parameters), signif(phi / omega,dig))


## ----nls-fitted-values, echo=FALSE--------------------------------------------
assignfun(nls.parameters, text=".nls")
tpeak.nls <- with(as.list(nls.parameters), signif(phi / omega,dig))
cc <- suppressMessages(confint(nlsfit))
## would like to use built-in knitr/Sexpr{} magic but ...
## from https://stackoverflow.com/questions/8366324/r-sweave-formatting-numbers-with-sexpr-in-scientific-notation
format_sn0 <- function(x, digits = dig,
                      max_ord = 2,
                      ensuremath = TRUE) {

    if (x==0) return("0")
    ord <- floor(log(abs(x), 10))
    if (abs(ord) <= abs(max_ord)) return(as.character(signif(x, dig)))
    x <- x / 10^ord
    x <- signif(x, digits=dig)
    r <- sprintf("%s\\times 10^{%s}", x, ord)
    if (ensuremath) r <- sprintf("\\ensuremath{%s}", r)
    return(r)
}
format_sn <- Vectorize(format_sn0,
                       vectorize.args = "x")
## format_sn <- knitr::knit_hooks$get("inline")
ci_fmt  <- function(p, cc, sn = TRUE, ...) {
    ffun <- if (sn) format_sn else format
    ( ffun(cc[p,], ...)
    |> paste(collapse = ", ")
    |> sprintf(fmt = "(%s)")
    )
}
assignfun(cc[,"2.5%"], text=".lwr")
assignfun(cc[,"97.5%"], text=".upr")


## ----invert-KM-analytically, echo=FALSE---------------------------------------
invert_params <- function(I0, params) {
    with(as.list(params), {
        tpeak <- phi / omega
        Reff <- 1 + 2 * omega * I0 * sinh(phi) * cosh(phi)/a
	q <- (Reff-1)/tanh(phi)
        gamma <- 2 * omega/q ## gamma has units of 1/week
	S0 <- 2 * Reff^2 * I0 * sinh(phi)^2/(Reff - 1)^2
	beta <- Reff * gamma/S0
	c(tpeak=tpeak, Reff=Reff, S0=S0, gamma=gamma, Tg=1/gamma, beta=beta)
    })
}
N.KM <- 10^6
I0.KM <- 1
orig.params.KM <- invert_params(I0 = I0.KM, params = KM.parameters)
orig.params.KM <- signif(orig.params.KM,dig)
assignfun(orig.params.KM)
R0.KM <- Reff.KM/S0.KM*N.KM
I0.nls <- 1
orig.params.nls <- invert_params(I0 = I0.nls, params = coef(nlsfit))
orig.params.nls <- signif(orig.params.nls,dig)
assignfun(orig.params.nls, ".nls")
R0.nls <- Reff.nls/S0.nls*N.KM


## ----delta-method-ci, echo=FALSE----------------------------------------------
convert_express <- list(
    tpeak=expression(phi/omega),
    Reff=expression((1 + 2 * omega * I0 * sinh(phi) * cosh(phi)/a)),
    gamma=expression(2 * omega/((1 + 2 * omega * I0 * sinh(phi) * cosh(phi)/a)-1)/tanh(phi)),
    S0=expression(2 * ((1 + 2 * omega * I0 * sinh(phi) * cosh(phi)/a))^2 * I0 * sinh(phi)^2/(((1 + 2 * omega * I0 * sinh(phi) * cosh(phi)/a)) - 1)^2),
    beta=expression((2 * omega*tanh(phi))*((1 + 2 * omega * I0 * sinh(phi) * cosh(phi)/a) - 1)/(2 * (1 + 2 * omega * I0 * sinh(phi) * cosh(phi)/a) * I0 * sinh(phi)^2)),
    Tg=expression(1/(2 * omega/((1 + 2 * omega * I0 * sinh(phi) * cosh(phi)/a)-1)/tanh(phi))),
    R0=expression(((1 + 2 * omega * I0 * sinh(phi) * cosh(phi)/a) - 1)^2/(2 * (1 + 2 * omega * I0 * sinh(phi) * cosh(phi)/a) * I0 * sinh(phi)^2)*N)
)

KMcifun <- function(fit,
                    I0,
                    N) {
    param <- coef(fit)
    ##dg/dx
    dgdx <- sapply(convert_express, function(x) {
        dd <- deriv(x, c("a", "omega", "phi"))
        ee <- eval(dd, c(as.list(param), I0=I0, N=N))
        attr(ee, "gradient")
    })

    vcov <- vcov(fit)
    est_vcov <- t(dgdx) %*% vcov %*% dgdx
    est_err <- sqrt(diag(est_vcov))
    est_par <- sapply(convert_express, eval, c(as.list(param), I0=I0, N=N))
    z <- -qnorm((1-0.95)/2)

    matrix(
        c(est_par - z * est_err, est_par + z * est_err),
        ncol=2,
        dimnames=list(
            c("tpeak", "Reff", "gamma", "S0", "beta", "Tg", "R0"),
            c("2.5%", "97.5%")
        )
    )
}

cc.delta <- KMcifun(nlsfit, I0=I0.KM, N=N.KM)
assignfun(cc.delta[,"2.5%"], text=".lwr")
assignfun(cc.delta[,"97.5%"], text=".upr")


## ----delta-method-band, echo=FALSE--------------------------------------------
traj_express <- expression(a * 1/cosh(omega*t - phi)^2)

traj_ci <- function(t,
                    parameters,
                    vcov,
                    level=0.95) {
    traj <- eval(deriv(traj_express, c("a", "omega", "phi")), c(as.list(parameters), t=t))

    traj_dgdx <- attr(traj, "gradient")

    est_vcov <- traj_dgdx %*% vcov %*% t(traj_dgdx)
    est_err <- sqrt(diag(est_vcov))
    ll <- (1-level)/2
    z <- -qnorm(ll)

    dd <- data.frame(
        times=t,
        estimate=traj,
        lwr=traj-est_err*z,
        upr=traj+est_err*z
    )
    names(dd)[3:4] <- c(paste(100*ll, "%"), paste(100*(1-ll), "%"))

    dd
}


## ----load fitode, message=FALSE-----------------------------------------------
library(fitode)


## ----sirmortmodel_correct-----------------------------------------------------
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


## ----sirmortfit_data, cache=TRUE, warning=FALSE-------------------------------
bombay2 <- rbind(
    c(times=bombay$week[1] -
          diff(bombay$week)[1], mort=NA),
    bombay
)


## ----sirmortfit,cache=TRUE,warning=FALSE,message=FALSE, tidy = FALSE----------
SIR_start <- c(beta=beta.nls, gamma=gamma.nls,
               I0=I0.KM, S0=S0.nls, k=50)
SIR_fit <- fitode( model = SIR_model, data = bombay2,
                  fixed = list(gamma=gamma.nls),
                  start = SIR_start, tcol = "week" )


## ----sirmortfitcoef, include=FALSE--------------------------------------------
coef.bombay.fitode.nb <- coef(SIR_fit)
ci.bombay.fitode.nb <- confint(SIR_fit)[,-1]
assignfun(coef.bombay.fitode.nb, ".bombay.fitode.nb")


## ----sirmortconfint, include=FALSE--------------------------------------------
ci.derived.bombay.fitode.nb <- confint(SIR_fit,
        parm=list(
            Reff~beta*S0/gamma.nls,
            R0~beta*N.KM/gamma.nls,
            Tg~1/gamma.nls
        ))
assignfun(ci.derived.bombay.fitode.nb[,1], ".bombay.fitode.nb")


## ----confband, include=FALSE--------------------------------------------------
tmax <- 33
tvals <- 0:tmax
SIR_confband <- predict(SIR_fit, times=tvals, level=0.95)$mort


## ----sirmortmodel, include=FALSE----------------------------------------------
SIR_model_ols <- odemodel(
    name="SIR model",
    model=list(
        S ~ - beta * S * I,
        I ~ beta * S * I - gamma * I,
        R ~ gamma * I
    ),
    observation = list(
        mort ~ ols(mean = R)
    ),
    initial=list(
        S ~ S0,
        I ~ I0,
        R ~ 0
    ),
    par=c("beta", "gamma", "S0", "I0"),
    diffnames="R"
)


## ----sirmortfit_correct, cache=TRUE, warning=FALSE, include=FALSE-------------
SIR_start_ols <- coef(SIR_fit)[-5]

SIR_fit_ols <- fitode(
    SIR_model_ols,
    data = bombay2,
    start = SIR_start_ols,
    fixed = list(gamma=gamma.nls),
    tcol = "week"
)


## ----sirmortfitcoefols, include=FALSE-----------------------------------------
coef.bombay.fitode.ols <- coef(SIR_fit_ols)
ci.bombay.fitode.ols <- confint(SIR_fit_ols)[,-1]
assignfun(coef.bombay.fitode.ols, ".bombay.fitode.ols")

ci.derived.bombay.fitode.ols <- confint(SIR_fit_ols,
        parm=list(
            Reff~beta*S0/gamma.nls,
            R0~beta*N.KM/gamma.nls,
            Tg~1/gamma.nls
        ))

assignfun(ci.derived.bombay.fitode.ols[,1], ".bombay.fitode.ols")


## ----sirmortfit_ols_confband, include=FALSE-----------------------------------
SIR_ols_confband <- predict(SIR_fit_ols, times=tvals, level=0.95)$mort


## ----bombay_likelihood_fit, cache=TRUE, warning=FALSE, message=FALSE, include=FALSE----
if (file.exists("bombay-likelihood.rda")) {
    load("bombay-likelihood.rda")
} else {
    nfits <- 100
    gammavec <- seq(1, 28, length.out=nfits)
    loglikvec <- rep(NA, nfits)
    fitlist <- vector('list', length(gammavec))

    for (i in 1:nfits) {
        if (i==1) {
            start_ll <- coef(SIR_fit2)
        } else {
            start_ll <- coef(fitlist[[which.max(loglikvec[1:(i-1)])]])
        }

        fitlist[[i]] <- try(fitode(
            SIR_model,
            data = bombay2,
            start = start_ll,
            fixed = list(gamma=gammavec[i]),
            tcol = "week"
        ))

        if (!inherits(fitlist[[i]], "try-error")) {
            loglikvec[i] <- logLik(fitlist[[i]])
        }
    }
    save("gammavec", "loglikvec", "fitlist", file="bombay-likelihood.rda")
}


## ----bombay_likelihood_plot, cache=TRUE, include=FALSE------------------------
plot(gammavec, loglikvec,
     xlab="gamma (1/week)",
     ylab="log likelihood",
     ylim=c(-140, -135.5))


## ----philapop, echo=TRUE, include=FALSE---------------------------------------
## https://en.wikipedia.org/wiki/Demographics_of_Philadelphia
philapop <- data.frame(year = c(1910,1920), population = c(1549008, 1823779))
##' @param pop data frame
##' @param year year to interpolate to
pop_interp <- function( pop, year ) {
    slope <- (pop[2,2] - pop[1,2])/(pop[2,1] - pop[1,1])
    out <- round( pop[1,2] + (year-pop[1,1])*slope )
    return(out)
}
philapop1918 <- pop_interp( philapop, 1918 )
print(philapop1918)
## reduce this by case fatality proportion so we have
## the initial susceptible population who will be recorded
## as a death if they get infected:
philacfp <- 0.025
philaS0 <- round(philapop1918 * philacfp)
print(philaS0)


## ----read.phila.data, include=FALSE-------------------------------------------
## we have to cut the data when we're using nbinom
## because there are long trailing zeroes..
## same with the beginning
phila1918a <- subset(phila1918, as.Date("1918-09-10") < date &
                                date < as.Date("1918-11-18"))
phila1918b <- data.frame(
    date=c(as.Date("1918-09-10"), phila1918a$date),
    mort=c(NA, phila1918a$mort)
)
## avoids dependence on lubridate
## not reliable for leap years!
mk_dec_date <- function(x) {
    p <- as.POSIXlt(x)
    1900 + p$year + p$yday/365
}
phila1918b$time <- mk_dec_date(phila1918b$date)
if (requireNamespace("lubridate")) {
    stopifnot(all.equal(
        phila1918b$time, lubridate::decimal_date(phila1918b$date)))
}



## ----phila.fit,cache=TRUE,warning=FALSE, include=FALSE------------------------
## time is in the unit of years
## so gamma=52 corresponds to a mean of 1 week
## not-so-random starting values chosen by hand
start.phila <- c(beta=4e-3, gamma=52*2,
                 S0=philaS0*0.9,
                 I0=3, k=50)
names(start.phila) <- c("beta","gamma","S0","I0", "k")
print(start.phila)

phila_fit <- fitode(
    SIR_model,
    data = phila1918b,
    start = start.phila,
    tcol = "time"
)
ppp <- predict(phila_fit, level=0.95)[[1]]


## ----philaassign, include=FALSE-----------------------------------------------
coef.phila.fitode <- coef(phila_fit)
ci.phila.fitode <- confint(phila_fit)[,-1]
assignfun(coef.phila.fitode, ".phila.fitode")


## ----sirphilaconfint, include=FALSE-------------------------------------------
ci.derived.phila.fitode <- confint(phila_fit,
        parm=list(
            Reff~beta*S0/gamma,
            R0~beta*philaS0/gamma,
            Tg~1/gamma
        ))
assignfun(ci.derived.phila.fitode[,1], ".phila.fitode")
philaparms <- as.list(c(coef.phila.fitode, ci.derived.phila.fitode[,1]))


## ----philaKMparms, include=FALSE----------------------------------------------
## get KM approx parms associated with fitode fit:
a.phila <- with(philaparms, a_fun(Reff=Reff, gamma=gamma, S0=S0, I0=I0))
omega.phila <- with(philaparms, omega_fun(Reff=Reff, gamma=gamma, S0=S0, I0=I0))
phi.phila <- with(philaparms, phi_fun(Reff=Reff, gamma=gamma, S0=S0, I0=I0))
KM.approx.phila <- c(a=a.phila, omega=omega.phila, phi=phi.phila)
print(KM.approx.phila)


## ----adjust philaKMparms, include=FALSE---------------------------------------
## adjust because above is giving nls problems:
KM.start.phila <- c(a=a.phila, omega=omega.phila, phi=phi.phila)
print(KM.start.phila)


## ----phila.nlsfit, include=FALSE----------------------------------------------
## shift time to start at 0 before passing to nls
df <- phila1918b[,c("time","mort")]
df[,"time"] <- df[,"time"] - df[1,"time"]
phila.nlsfit <- nls(mort ~ KM_approx(time, a, omega, phi),
              data = df,
              start = KM.start.phila)
phila.nls.parameters <- coef(phila.nlsfit)
ci.phila.nls <- confint(phila.nlsfit)
assignfun(phila.nls.parameters, ".phila.nls")
print(phila.nls.parameters)


## ----philanlsR0, include=FALSE------------------------------------------------
cc.phila.delta <- KMcifun(phila.nlsfit, I0=I0.phila.fitode, N=philaS0)
orig.params.phila <- invert_params(I0=I0.phila.fitode, coef(phila.nlsfit))
orig.params.phila <- signif(orig.params.phila,dig)
assignfun(orig.params.phila, text = ".phila.nls")
R0.phila.nls <- Reff.phila.nls/S0.phila.nls*philaS0


## ----compute-stochastic-SIR-via-sirr, include=FALSE---------------------------
sirstoch_fn <- system.file("vignette_data", "sirstoch.RData", package = "fitode")
if (file.exists(sirstoch_fn)) {
    load(sirstoch_fn)
} else {
    library(sirr)
    R0 <- 5
    N <- 2000
    mm <- create_SIRmodel(R0=R0, N=N)
    ii <- set_inits(mm)
    iii <- ii
    iii[2] <- exp(iii[2])
    names(iii) <- c("S","I","R")
    ## FIX: this should not be necessary, but it is a bug in sirr
    ##      associated with N != 1.
    mm1 <- create_SIRmodel(R0=R0, N=1)
    sir.det <- compute_SIRts(mm1)
    ## get associated params for KM approx:
    a.det <- a_fun(Reff=R0*iii["S"]/2000, gamma=mm$gamma, S0=iii["S"], I0=iii["I"])
    omega.det <- omega_fun(Reff=R0*iii["S"]/2000, gamma=mm$gamma, S0=iii["S"], I0=iii["I"])
    phi.det <- phi_fun(Reff=R0*iii["S"]/2000, gamma=mm$gamma, S0=iii["S"], I0=iii["I"])
    ##
    tt <- 0:10
    sir.stoch <- get_ssa_soln(mm, stopcrit = NULL, inits=ii, times=tt)
    save(R0, N, mm, ii, iii, mm1, sir.det, a.det, omega.det, phi.det, tt,
         sir.stoch, file="sirstoch.RData")
}


## ----colours, include=FALSE---------------------------------------------------
my.red <- "#E41A1C"     # colour-blind friendly red
my.blue <- "#377EB8"    # colour-blind friendly blue
my.green <- "#4DAF4A" # colour-blind friendly green
my.yellow <- "#FFD92F" # colour-blind friendly yellow
my.orange <- "#D55E00" # colour-blind friendly orange


## ----sirstochdata, include=FALSE----------------------------------------------
## convert to daily mortality
mort.stoch <- diff(sir.stoch$R[!duplicated(ceiling(sir.stoch$time*7),fromLast = TRUE)])
mort.time <- ceiling(sir.stoch$time[!duplicated(ceiling(sir.stoch$time*7),fromLast = TRUE)]*7)
mort.which <- which(diff(mort.time)==1)

sir.stoch.data <- data.frame(
    time=mort.time[mort.which]/7,
    mort=c(NA, mort.stoch[head(mort.which, -1)])
)


## ----sirnlsfit, cache=TRUE, warning=FALSE, include=FALSE----------------------
nlsfit.stoch <- nls(mort ~ KM_approx(time, a, omega, phi),
              data = sir.stoch.data,
              start = c(a = unname(a.det), omega = unname(omega.det),
                        phi = unname(phi.det)))

nls.stoch.parameters <- coef(nlsfit.stoch)
assignfun(nls.stoch.parameters, text=".nls.stoch")

orig.params.KM.stoch <- invert_params(I0 = exp(ii[[2]]), params = c(a=unname(a.det), phi=unname(phi.det), omega=unname(omega.det)))
orig.params.KM.stoch <- signif(orig.params.KM.stoch)
assignfun(orig.params.KM.stoch, ".KM.stoch")

orig.params.nls.stoch <- invert_params(I0 = exp(ii[[2]]), params = coef(nlsfit.stoch))
orig.params.nls.stoch <- signif(orig.params.nls.stoch)
assignfun(orig.params.nls.stoch, ".nls.stoch")

cc.nls.stoch <- suppressMessages(confint(nlsfit.stoch))

cc.delta.stoch <- KMcifun(nlsfit.stoch, I0 = exp(ii[[2]]), N=N)
assignfun(cc.delta.stoch[,"2.5%"], text=".lwr.nls.stoch")
assignfun(cc.delta.stoch[,"97.5%"], text=".upr.nls.stoch")



## ----sirstochfit, cache=TRUE, warning=FALSE, include=FALSE--------------------
start.stoch <- c(beta=mm$beta, gamma=mm$gamma,
                 S0=iii["S"], I0=iii["I"])
names(start.stoch) <- c("beta","gamma","S0","I0")
print(start.stoch)

SIR_stoch_fit <- fitode(
    SIR_model,
    data = sir.stoch.data,
    start = c(start.stoch, k=10),
    # fixed = c(gamma = 1), # Tg = 1/gamma = 1 week
    tcol = "time"
)
sss <- predict(SIR_stoch_fit, level=0.95)[[1]]
summary(SIR_stoch_fit)


## ----stochestassign, include=FALSE--------------------------------------------
coef.stoch.fitode <- coef(SIR_stoch_fit)
ci.stoch.fitode <- confint(SIR_stoch_fit)[,-1]
assignfun(coef.stoch.fitode, ".stoch.fitode")


## ----sirstochconfint, include=FALSE-------------------------------------------
ci.derived.stoch.fitode <- confint(SIR_stoch_fit,
        parm=list(
            Reff~beta*S0/gamma,
            R0~beta*2000/gamma
        ))
assignfun(ci.derived.stoch.fitode[,1], ".stoch.fitode")


## ----figfuns, echo=FALSE------------------------------------------------------
lwd <- 5
use_OI <- TRUE
pp <- palette.colors()
if (FALSE) plot(1:8, 1:8, cex = 5, pch = 16, col = pp)
if (use_OI) {
    ## these were chosen from the O-I palette to approximately
    ## match the previous choices (substituted an orange for yellow,
    ## as I found the yellow a bit too light)
    col.KM <- pp[6]
    col.nls <- pp[7]
    col.fitode.nb <- pp[2]
    col.fitode.ols <- pp[3]
    col.determ <- pp[4]
} else {
    col.KM <- my.blue
    col.nls <- my.red
    col.fitode.nb <- my.yellow
    col.fitode.ols <- "grey80" # "orange" # my.orange
    col.determ <- my.green
}


transparent_colour <- function(col,alpha=150/255) {
    alpha <- round(alpha * 255)
    v <- col2rgb(col)[,1] # color as rgb vector
    tcol <- rgb(v["red"],v["green"],v["blue"],alpha=alpha,maxColorValue = 255)
    return(tcol)
}
col.fitodenbCI <- transparent_colour(col.fitode.nb, alpha=0.4)
col.fitodeCI <- transparent_colour(col.fitode.ols, alpha=0.4)
col.nlsCI <- transparent_colour(col.nls, alpha=0.4)
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
draw_legend <- function(lwd=5, pt.bg="white") {
    legend("topleft", bty="n", lwd=c(2,lwd,lwd,lwd*0.60),
       lty=c(NA,"solid","solid","solid"),
       col=c("black",col.KM,col.nls,col.fitode.nb),
       pch=c(21,NA,NA,NA),
       pt.bg=c(pt.bg,NA,NA,NA),
       legend=c("observed data","KM","nls","fitode (nbinom)"))
}


## ----testing, include=FALSE, results = "hide"---------------------------------
format_sn0(cc["a", 1])
ci_fmt("a", cc)


## ----Bombay-figure, echo=FALSE, fig.height=5, fig.show="hold", out.width="95%"----
tvals <- seq(0, tmax, length.out=1000)
setup_plot(ylim = c(0, 1000), xlim=c(0,tmax),
         xlab="weeks", ylab="plague deaths")
draw_confband(traj_ci(tvals, nls.parameters, vcov(nlsfit)), col=col.nlsCI)
##DE: skipping fitode confband because there is too much overlap:
##draw_confband(SIR_confband, col=col.fitodenbCI)
curve(KM_approx(x, a=a.KM, omega=omega.KM, phi=phi.KM), from=0, to=tmax, n=1000, add=TRUE,
      lwd=lwd, xpd=NA, col=col.KM)
curve(KM_approx(x, a=a.nls, omega=omega.nls, phi=phi.nls), from=0, to=tmax, n=1000, add=TRUE,
      lwd=lwd, xpd=NA, col=col.nls)
lines(estimate ~ times, data=SIR_confband, col=col.fitode.nb, lwd=lwd*0.60, lty="dotted")
points(mort ~ week, data = fitode::bombay, xpd=NA, pch=21, bg="white", lwd=2)
legend("topleft", bty="n", lwd=c(2,lwd,lwd,lwd*0.60),
       lty=c(NA,"solid","solid","dotted"),
       col=c("black",col.KM,col.nls,col.fitode.nb),
       pch=c(21,NA,NA,NA),
       pt.bg=c("white",NA,NA,NA),
       legend=c("observed data","KM's (1927) fit","our nls fit","fitode (nbinom)"))

setup_plot(ylim=c(0, 1000), xlim=c(0,tmax), #DE: ymin was 5, not sure why
         xlab="weeks", ylab="plague deaths"
         ##, at=c(25, 50, 100, 200, 400, 800)
         )
draw_confband(SIR_confband, col=col.fitodenbCI)
draw_confband(SIR_ols_confband, col=col.fitodeCI)
lines(estimate ~ times, data=SIR_confband, col=col.fitode.nb, lwd=lwd)
lines(estimate ~ times, data=SIR_ols_confband, col=col.fitode.ols, lwd=lwd*0.60)
curve(KM_approx(x, a=a.nls, omega=omega.nls, phi=phi.nls), from=0, to=tmax, n=250, add=TRUE,
      lwd=lwd*0.60, xpd=NA, col=col.nls, lty="dotted")
points(mort ~ week, data = fitode::bombay, xpd=NA, pch=21, bg="white", lwd=2)
legend("topleft", bty="n", lwd=c(2,lwd,lwd*0.60,lwd*0.6),
       lty=c(NA,"solid","solid","dotted"),
       col=c("black",col.fitode.nb,col.fitode.ols,col.nls),
       pch=c(21,NA,NA,NA),
       pt.bg=c("white",NA,NA,NA),
       legend=c("observed data","fitode (nbinom)", "fitode (ols)", "our nls fit"))


## ----phila-figure, echo=FALSE, fig.height=5-----------------------------------
setup_plot(ylim = c(0, 800), xlim = range(phila1918b$time), ylab="daily P&I deaths",
           xaxt='n')
axis(side=1, at=phila1918b$time[c(1, 31, 62)],
     labels=c("Sep", "Oct", "Nov"))
draw_confband(ppp)
tvals <- with(phila1918b, seq(min(time),max(time), length=1000))
t0 <- min(tvals)
tt_phila <- traj_ci(tvals-t0, phila.nls.parameters, vcov(phila.nlsfit))

tt_phila$times <- tt_phila$times + t0

draw_confband(tt_phila, col=col.nlsCI)
curve(KM_approx(x-t0, a=a.phila, omega=omega.phila, phi=phi.phila)/365, from=t0, to=max(tvals), n=1000, add=TRUE, xpd=NA, col=col.fitode.nb, lty="dotted", lwd=lwd*0.60)
lines(tvals,
      with(as.list(phila.nls.parameters),KM_approx(tvals-t0, a=a, omega=omega, phi=phi)),
      xpd=NA, lwd=lwd, col=col.nls)
lines(estimate ~ times, data=ppp, col=col.fitode.nb, lwd=lwd*0.8, lty="solid")
points(mort ~ time,
       data = phila1918b, xpd=NA, pch=21, bg="white", lwd=2)
legend("topleft", bty="n", lwd=c(2,lwd,lwd*0.8,lwd*0.60),
       lty=c(NA,"solid","solid","dotted"),
       col=c("black",col.nls,col.fitode.nb,col.fitode.nb),
       pch=c(21,NA,NA,NA),
       pt.bg=c("white",NA,NA,NA),
       legend=c("observed data","nls fit of KM approx",
                "fitode (nbinom)","KM approx to fitode fit"))


## ----plot_stochastic_SIR, echo=FALSE, fig.height=5, fig.show="hold", out.width="95%"----
my.tmax <- 8
plot(NA, NA, bty="L", xlab="", ylab="", las=1, xaxs="i", yaxs="i",
     ##xlim=range(tt),
     xlim=c(0,my.tmax), # FIX: hack: max(tt) in sirstoch.RData is 10, which is too large
     ylim=range(sir.stoch.data$mort, na.rm=TRUE))
draw_confband(sss, col=col.fitodenbCI)
lines(estimate ~ times, data=sss, col=col.fitode.nb, lwd=lwd)
lines(sir.det$tau, 1 * sir.det$I * N/7, col=col.determ, lwd=lwd*0.8) ## dR/dt
curve(KM_approx(x, a=a.det, omega=omega.det, phi=phi.det)/7, from=0, to=my.tmax, n=1000, add=TRUE,
      xpd=NA, col=col.determ, lty="dotted", lwd=lwd*0.6)
points(mort ~ time, data = sir.stoch.data, xpd=NA, pch=21, bg="white", lwd=2)
legend("topright", bty="n", lwd=c(2,lwd,lwd*0.8,lwd*0.6),
       lty=c(NA,"solid","solid","dotted"),
       col=c("black",col.fitode.nb,col.determ,col.determ),
       pch=c(21,NA,NA,NA),
       pt.bg=c("white",NA,NA,NA),
       ##adj = c(0, 0.8),  ## shift text so first line aligns with guide line
       legend=c("simulated data","fitode (nbinom)","deterministic","KM approx"))

plot(NA, NA, bty="L", xlab="", ylab="", las=1, xaxs="i", yaxs="i",
     ##xlim=range(tt),
     xlim=c(0,my.tmax), # FIX: hack: max(tt) in sirstoch.RData is 10, which is too large
     ylim=range(sir.stoch.data$mort, na.rm=TRUE))
draw_confband(traj_ci((0:(my.tmax*7))/7, nls.stoch.parameters, vcov(nlsfit.stoch)), col=col.nlsCI)
draw_confband(sss, col=col.fitodenbCI)
lines(estimate ~ times, data=sss, col=col.fitode.nb, lwd=lwd)
curve(KM_approx(x, a=a.nls.stoch, omega=omega.nls.stoch, phi=phi.nls.stoch),
      from=0, to=my.tmax, n=1000, add=TRUE,
      xpd=NA, col=col.nls, lwd=lwd*0.8, lty="solid")
points(mort ~ time, data = sir.stoch.data, xpd=NA, pch=21, bg="white", lwd=2)
legend("topright", bty="n", lwd=c(2,lwd,lwd*0.8),
       lty=c(NA,"solid","solid"),
       col=c("black",col.fitode.nb,col.nls),
       pch=c(21,NA,NA),
       pt.bg=c("white",NA,NA),
       legend=c("simulated data","fitode (nbinom)","nls fit of KM approx"))

