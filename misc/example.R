library(fitode)

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
    par=c("beta", "gamma", "N", "i0", "size1", "size2"),
    link=c(i0="logit")
)

parms <-c(beta=1,gamma=0.5, N=1000, i0=1e-2, size1=10, size2=10)

df <- simulate(SI_model, times=1:20, parms=parms, seed=101)[,c("time", "susceptible", "infected")]

system.time(ff <- fitode(
    model=SI_model,
    data=df,
    start=parms,
    tcol="time"
))

plot(ff, level=0.95)

confint(ff, parms=c("N", "gamma"))
confint(ff, parms=c("N", "gamma"), method="profile")
confint(ff, parms=c("N", "gamma"), method="wmvrnorm")

SI_model_c <- new("model.ode",
    name = "SI",
    model = list(
        S ~ - beta*S*I/(S+I),
        I ~ beta*S*I/(S+I) - gamma*I,
        R ~ gamma*I
    ),
    observation <- list(
        cases ~ dnbinom(mu=R, size=size)
    ),
    initial = list(
        S ~ N * s0,
        I ~ N * i0,
        R ~ 0
    ),
    par=c("beta", "gamma", "N", "i0", "s0", "size"),
    diffnames=c("R")
)

harbin_1910a <- rbind(
    c(0, NA),
    harbin_1910
)

start <- c(beta=0.6, gamma=0.5, i0=2.5e-4, s0=0.13, size=5)

system.time(ff4 <- fitode(
    model=SI_model_c,
    start=start,
    data=harbin_1910a,
    link = list(
        i0="logit",
        s0="logit"
    ),
    fixed=c(N=25000)
))

plot(ff4, level=0.95)
plot(ff4, method="wmvrnorm", level=0.95)
