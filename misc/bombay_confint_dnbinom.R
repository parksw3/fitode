library(fitode)
library(dplyr)
library(ggplot2); theme_set(theme_bw(base_family="Times"))
library(egg)
library(gridExtra)

if (file.exists("bombay_confint_dnbinom.rda")) {
    load("bombay_confint_dnbinom.rda")
    nfits <- 500
    gammavec <- seq(1, 100, length.out=nfits)
} else {
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
    save("fitlist_gamma", file="bombay_confint_dnbinom.rda")
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
    geom_point(aes(gamma, logLik - max(logLik))) +
    scale_x_continuous(expression(gamma)) +
    scale_y_continuous("Profile log-likelihood") +
    ggtitle("A")

g2 <- ggplot(gammadata) +
    geom_point(aes(gamma, beta*1e5)) +
    labs(x = expression("recovery rate"~(gamma)),
         y = expression("transmission rate"~(beta)~(""%*%10^5))) + 
## scale_x_continuous(expression(gamma)) +
##    scale_y_continuous(expression(beta)) +
    ggtitle("B")

g3 <- ggplot(gammadata) +
    geom_point(aes(gamma, S0/1e5)) +
    labs(x = expression("recovery rate"~(gamma)),
         y = expression("init susceptibles"~S(0)~(""/10^5)))
    ggtitle("C")

g4 <- ggplot(gammadata) +
    geom_point(aes(gamma, I0))  +
    labs(x = expression("recovery rate"~(gamma)),
         y = expression("init infectives"~I(0))) +
    ggtitle("D")

g5 <- ggplot(gammadata) +
    geom_point(aes(gamma, k)) +
    labs(x = expression("recovery rate"~(gamma)),
         y = expression("dispersion"~(k))) +
    ggtitle("E")
gcomb <- ggarrange(g2, g3, g4, g5)

gfinal <- arrangeGrob(g1, gcomb, nrow=1, widths=c(1, 2))

ggsave("bombay_confint_dnbinom.pdf", gfinal, width=8, height=4)
