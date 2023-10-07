library(fitode)
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

SIR_start <- c(beta = 7.83e-05, gamma = 4.11, I0 = 1, S0 = 57400, k = 50)
SIR_fit <- fitode(
    model = SIR_model,
    data = bombay2,
    start = SIR_start,
    fixed = list(gamma=SIR_start[["gamma"]]),
    tcol = "week"
)

nfits <- 20
betavec <- seq(7e-5, 11e-5, length.out=nfits)
S0vec <- seq(4.5e4, 5.5e4, length.out=nfits)
I0vec <- seq(0.5, 1.2, length.out=nfits)
kvec <- seq(20, 120, length.out=nfits)

loglikvec_beta <- loglikvec_S0 <- loglikvec_I0 <- loglikvec_k <- rep(NA, nfits)
fitlist_beta <- fitlist_S0 <- fitlist_I0 <- fitlist_k <- vector('list', nfits)

for (i in 1:nfits) {
    print(i)
    if (i==1) {
        start_ll <- coef(SIR_fit)
    } else {
        start_ll <- coef(fitlist_beta[[which.max(loglikvec_beta[1:(i-1)])]])
    }

    ff <- try(fitode(
        SIR_model,
        data = bombay2,
        start = start_ll,
        fixed = list(gamma=SIR_start[["gamma"]], beta=betavec[[i]]),
        tcol = "week",
        skip.hessian = TRUE
    ))

    if (!inherits(ff, "try-error")) {
        fitlist_beta[[i]] <- ff
        loglikvec_beta[i] <- logLik(ff)
    }
}

for (i in 1:nfits) {
    print(i)
    if (i==1) {
        start_ll <- coef(SIR_fit)
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

for (i in 1:nfits) {
    print(i)
    if (i==1) {
        start_ll <- coef(SIR_fit)
    } else {
        start_ll <- coef(fitlist_k[[which.max(loglikvec_k[1:(i-1)])]])
    }

    ff <- try(fitode(
        SIR_model,
        data = bombay2,
        start = start_ll,
        fixed = list(gamma=SIR_start[["gamma"]], k=kvec[[i]]),
        tcol = "week",
        skip.hessian = TRUE
    ))

    if (!inherits(ff, "try-error")) {
        fitlist_k[[i]] <- ff
        loglikvec_k[i] <- logLik(ff)
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

data_k <- data.frame(
    k=kvec,
    logLik=loglikvec_k
)

g1 <- ggplot(data_beta) +
    geom_point(aes(beta, logLik)) +
    scale_x_continuous(limits=c(8e-5, 10.2e-5)) +
    scale_y_continuous("Profile log-likelihood", limits=c(-144, -137)) +
    theme(
        panel.grid = element_blank()
    )

g2 <- ggplot(data_S0) +
    geom_point(aes(S0, logLik)) +
    scale_x_continuous(limits=c(NA, 53000)) +
    scale_y_continuous("Profile log-likelihood", limits=c(-144, -137)) +
    theme(
        panel.grid = element_blank()
    )

g3 <- ggplot(data_I0) +
    geom_point(aes(I0, logLik)) +
    scale_y_continuous("Profile log-likelihood", limits=c(-144, -137)) +
    theme(
        panel.grid = element_blank()
    )

g4 <- ggplot(data_k) +
    geom_point(aes(k, logLik)) +
    scale_x_continuous(limits=c(NA, 126)) +
    scale_y_continuous("Profile log-likelihood", limits=c(-144, -137)) +
    theme(
        panel.grid = element_blank()
    )

gcomb <- ggarrange(g1, g2, g3, g4, draw=FALSE)

save("fitlist_beta", "fitlist_S0", "fitlist_I0", "fitlist_k", file="bombay_confint_dnbinom_fixgamma.rda")
ggsave("bombay_confint_dnbinom_fixgamma.pdf", gcomb, width=8, height=6)
