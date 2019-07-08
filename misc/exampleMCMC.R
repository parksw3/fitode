library(fitode)

model <- new("model.ode",
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

## weird prior just for testing
prior <- list(
    beta~dgamma(shape=2, rate=1)
)

parms <-c(beta=1,gamma=0.5, N=1000, i0=1e-2, size1=10, size2=10)

df <- simulate(model, times=1:20, parms=parms, seed=101)[,c("times", "susceptible", "infected")]

system.time(ff <- fitodeMCMC(
    model=model,
    data=df,
    start=parms,
    vcov = 1*diag(length(parms)),
    prior=prior,
    chains=2, thin=10,
    prior.only=TRUE
))

plot(as.mcmc.list(ff[[1]])[,1])

