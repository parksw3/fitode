library(fitode)

harbin_1910 <- data.frame(
    times=1:68,
    cases=c(8,3,6,10,9,8,5,12,11,
            13,14,13,9,10,15,22,
            17,15,20,31,22,50,29,
            43,29,42,43,82,58,22,
            42,73,78,63,54,69,83,
            100,99,70,100,98,100,
            85,96,98,96,105,64,65,
            75,74,86,74,76,62,82,
            60,68,55,36,34,38,28,
            23,16,13,14)
)

SI_model <- new("model.ode",
    name = "SI",
    model = list(
        S ~ - beta*S*I/N,
        I ~ beta*S*I/N - gamma*I
    ),
    observation = list(
        susceptible ~ dnbinom(mu=S, size=size1),
        infected ~ dnbinom(mu=I, size=size2)
    ),
    initial = list(
        S ~ N * (1 - i0),
        I ~ N * i0
    ),
    par=c("beta", "gamma", "N", "i0", "size1", "size2")
)

tvec <- 1:20

f <- ode.solve(SI_model, tvec, c(beta=1,gamma=0.5, N=1000, i0=1e-2, sigma1=0.1, sigam2=0.1))

set.seed(101)
df <- data.frame(
    day=tvec,
    susceptible=rnbinom(length(tvec), mu=f@solution$S, size=10),
    infected=rnbinom(length(tvec), mu=f@solution$I, size=10)
)

start <- c(beta=1, gamma=0.5, N=1000, i0=1e-2, size1=10, size2=10)

system.time(ff <- fitode(
    model=SI_model,
    data=df,
    start=start,
    tcol="day",
    link = list(
        i0="logit"
    )
))

SI_model_c <- new("model.ode",
    name = "SI",
    model = list(
        S ~ - beta*S*I/N,
        I ~ beta*S*I/N - gamma*I,
        R ~ gamma*I
    ),
    observation <- list(
        cases ~ dpois(lambda=R)
    ),
    initial = list(
        S ~ N * (1 - i0),
        I ~ N * i0,
        R ~ 0
    ),
    par=c("beta", "gamma", "N", "i0"),
    diffnames=c("R")
)

harbin_1910a <- rbind(
    c(0, NA),
    harbin_1910
)

start <- c(beta=0.429, gamma=0.3, N=6340, i0=1.07e-3)

system.time(ff4 <- fitode(
    model=SI_model_c,
    start=start,
    data=harbin_1910a, tcol="day",
    link = list(
        i0="logit"
    )
))

plot(harbin_1910a$cases)
lines((ode.solve(SI_model_c, harbin_1910a$times, coef(ff4))@solution$R))


