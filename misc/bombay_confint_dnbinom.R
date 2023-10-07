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
        mort ~ dnbinom(mu = R+1e-301, size = k)
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

nfits <- 500
gammavec <- seq(1, 100, length.out=nfits)
loglikvec_gamma <- rep(NA, nfits)
fitlist_gamma <- vector('list', length(gammavec))

for (i in 1:nfits) {
    print(i)
    if (i==1) {
        start_ll <- coef(SIR_fit)
    } else {
        start_ll <- coef(fitlist_gamma[[which.max(loglikvec_gamma[1:(i-1)])]])
    }

    fitlist_gamma[[i]] <- try(fitode(
        SIR_model,
        data = bombay2,
        start = start_ll,
        fixed = list(gamma=gammavec[i]),
        tcol = "week",
        skip.hessian = TRUE
    ))

    if (!inherits(fitlist_gamma[[i]], "try-error")) {
        loglikvec_gamma[i] <- logLik(fitlist_gamma[[i]])
    }
}

gammadata <- lapply(fitlist_gamma, function(x) {
    as.data.frame(t(c(coef(x), logLik=logLik(x))))
}) %>%
    bind_rows %>%
    mutate(
        gamma=gammavec
    ) %>%
    filter(logLik > -140)

g1 <- ggplot(gammadata) +
    geom_point(aes(gamma, logLik)) +
    scale_y_continuous("Profile log-likelihood")

g2 <- ggplot(gammadata) +
    geom_point(aes(gamma, beta))

g3 <- ggplot(gammadata) +
    geom_point(aes(gamma, S0))

g4 <- ggplot(gammadata) +
    geom_point(aes(gamma, I0))

g5 <- ggplot(gammadata) +
    geom_point(aes(gamma, k))

gcomb <- ggarrange(g1, g2, g3, g4, g5, nrow=2)

save("fitlist_gamma", file="bombay_confint_dnbinom.rda")
ggsave("bombay_confint_dnbinom.pdf", gcomb, width=8, height=6)
