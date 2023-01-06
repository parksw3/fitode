bombay2 <- rbind(
    c(times=bombay$week[1] -
          diff(bombay$week)[1], mort=NA),
    bombay
)

nfits <- 21
gammavec <- seq(1, 14, length.out=nfits)
loglikvec <- rep(NA, nfits)
fitlist <- vector('list', length(gammavec))

for (i in 1:nfits) {
    if (i==1) {
        start_ll <- coef(SIR_fit2)
    } else {
        start_ll <- coef(fitlist[[i-1]])
    }

    ## hini=1/4 to solve the ode at a finer scale.. so that things don't go
    ## crazy when we're starting with very high gamma
    fitlist[[i]] <- fitode(
        SIR_model,
        data = bombay2,
        start = start_ll,
        fixed = list(gamma=gammavec[i]),
        tcol = "week"
    )

    loglikvec[i] <- logLik(fitlist[[i]])
}
