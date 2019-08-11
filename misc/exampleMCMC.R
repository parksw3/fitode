library(fitode)

model <- odemodel(
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

parms <- c(beta=1,gamma=0.5, N=1000, i0=1e-2, size1=10, size2=10)

df <- simulate(model, times=1:20, parms=parms, seed=101)

system.time(ff <- fitodeMCMC(
    model=model,
    data=df,
    start=parms,
    proposal.vcov = 0.01*diag(length(parms)),
    thin=10,
    chains=2,
    prior=list(
        beta~dgamma(shape=2, rate=2),
        gamma~dgamma(shape=2, rate=4),
        N~dgamma(shape=2, rate=1/500),
        i0~dbeta(shape1=1, shape2=1),
        size1~dgamma(shape=2, rate=1/5),
        size2~dgamma(shape=2, rate=1/5)
    )
))

plot(ff, level=0.95)

summary(ff)

lattice::xyplot(ff@mcmc)


coda::gelman.diag(ff@mcmc)

plot(ff@mcmc[,5])
