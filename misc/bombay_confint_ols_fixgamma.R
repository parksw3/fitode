library(fitode)
library(dplyr)
library(ggplot2); theme_set(theme_bw())
library(egg)

SIR_model <- odemodel(
    name="SIR model",
    model=list(
        S ~ - beta * S * I,
        I ~ beta * S * I - gamma * I,
        R ~ gamma * I
    ),
    observation = list(
        mort ~ dnorm(mean = R, sd=sigma)
    ),
    diffnames="R",
    initial=list(
        S ~ S0,
        I ~ I0,
        R ~ 0
    ),
    par=c("beta", "gamma", "S0", "I0", "sigma")
)

bombay2 <- rbind(
    c(times=bombay$week[1] -
          diff(bombay$week)[1], mort=NA),
    bombay
)

SIR_start <- c(beta = 7.83e-05, gamma = 4.11, I0 = 1, S0 = 57400, sigma=50)
SIR_fit <- fitode(
    model = SIR_model,
    data = bombay2,
    start = SIR_start,
    fixed = list(gamma=SIR_start[["gamma"]]),
    tcol = "week"
)

nfits <- 100
betavec <- seq(7e-5, 11e-5, length.out=nfits)
S0vec <- seq(4.5e4, 7.5e4, length.out=nfits)
I0vec <- seq(0.5, 2, length.out=nfits)
kvec <- seq(20, 120, length.out=nfits)

loglikvec_beta <- loglikvec_S0 <- loglikvec_I0 <- loglikvec_k <- rep(NA, nfits)
fitlist_beta <- fitlist_S0 <- fitlist_I0 <- fitlist_k <- vector('list', nfits)

for (i in 1:nfits) {
    print(i)
    if (i==1) {
        start_ll <- c(gamma=3.11, I0=1, S0=45000, sigma=70)
    } else {
        start_ll <- coef(fitlist_beta[[which.max(loglikvec_beta[1:(i-1)])]])
    }

    if (i==1) {
        ff <- try(fitode(
            SIR_model,
            data = bombay2,
            start = start_ll,
            fixed = list(gamma=SIR_start[["gamma"]], beta=betavec[[i]]),
            tcol = "week",
            skip.hessian = TRUE
        ))
    } else {
        ff <- try(fitode(
            SIR_model,
            data = bombay2,
            start = start_ll,
            fixed = list(gamma=SIR_start[["gamma"]], beta=betavec[[i]]),
            tcol = "week",
            skip.hessian = TRUE,
            method="BFGS"
        ))
    }


    if (!inherits(ff, "try-error")) {
        fitlist_beta[[i]] <- ff
        loglikvec_beta[i] <- logLik(ff)
    }
}

for (i in 1:nfits) {
    print(i)
    if (i==1) {
        start_ll <- c(beta=7e-5, gamma=3.11, I0=1, sigma=70)
    } else {
        start_ll <- coef(fitlist_S0[[which.max(loglikvec_S0[1:(i-1)])]])
    }

    ff <- try(fitode(
        SIR_model,
        data = bombay2,
        start = start_ll,
        fixed = list(gamma=SIR_start[["gamma"]], S0=S0vec[[i]]),
        tcol = "week",
        skip.hessian = TRUE
    ))

    if (!inherits(ff, "try-error")) {
        fitlist_S0[[i]] <- ff
        loglikvec_S0[i] <- logLik(ff)
    }
}

for (i in 1:nfits) {
    print(i)
    if (i==1) {
        start_ll <- coef(SIR_fit)
    } else {
        start_ll <- coef(fitlist_I0[[which.max(loglikvec_I0[1:(i-1)])]])
    }

    ff <- try(fitode(
        SIR_model,
        data = bombay2,
        start = start_ll,
        fixed = list(gamma=SIR_start[["gamma"]], I0=I0vec[[i]]),
        tcol = "week",
        skip.hessian = TRUE
    ))

    if (!inherits(ff, "try-error")) {
        fitlist_I0[[i]] <- ff
        loglikvec_I0[i] <- logLik(ff)
    }
}

data_beta <- data.frame(
    beta=betavec,
    logLik=loglikvec_beta
)

data_S0 <- data.frame(
    S0=S0vec,
    logLik=loglikvec_S0
)

data_I0 <- data.frame(
    I0=I0vec,
    logLik=loglikvec_I0
)

g1 <- ggplot(filter(data_beta, logLik > -170)) +
    geom_point(aes(beta, logLik)) +
    geom_smooth(aes(beta, logLik), se=FALSE, col="black", method="gam", lwd=0.7) +
    # scale_x_continuous(limits=c(8e-5, 10.2e-5)) +
    scale_y_continuous("Profile log-likelihood", limits=c(-170, -164)) +
    theme(
        panel.grid = element_blank()
    )

g2 <- ggplot(filter(data_S0, logLik > -170)) +
    geom_point(aes(S0, logLik)) +
    geom_smooth(aes(S0, logLik), se=FALSE, col="black", method="gam", lwd=0.7) +
    # scale_x_continuous(limits=c(NA, 53000)) +
    scale_y_continuous("Profile log-likelihood", limits=c(-170, -164)) +
    theme(
        panel.grid = element_blank()
    )

g3 <- ggplot(data_I0) +
    geom_point(aes(I0, logLik)) +
    geom_smooth(aes(I0, logLik), se=FALSE, col="black", method="gam", lwd=0.7) +
    scale_y_continuous("Profile log-likelihood", limits=c(-170, -164)) +
    theme(
        panel.grid = element_blank()
    )

gcomb <- ggarrange(g1, g2, g3, nrow=2, draw=FALSE)

save("fitlist_beta", "fitlist_S0", "fitlist_I0", file="bombay_confint_ols_fixgamma.rda")
ggsave("bombay_confint_ols_fixgamma.pdf", gcomb, width=8, height=6)
