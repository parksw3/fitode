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
## parameters are fitted on the log scale by default
model@link

parms <- c(beta=1,gamma=0.5, N=1000, i0=1e-2, size1=10, size2=10)

df <- simulate(model, times=1:20, parms=parms, seed=101)
matplot(df[,1],df[,-1],type="l",lty=1)

prior_list <- list(
    beta~dgamma(shape=2, rate=2),
    gamma~dgamma(shape=2, rate=4),
    N~dgamma(shape=2, rate=1/500),
    i0~dbeta(shape1=1, shape2=1),
    size1~dgamma(shape=2, rate=1/5),
    size2~dgamma(shape=2, rate=1/5)
)

## priors with mcmc
## i.e. get full Bayesian posterior distribution
system.time(ff <- fitodeMCMC(
    model=model,
    data=df,
    start=parms,
    proposal.vcov = 0.01*diag(length(parms)),
    thin=10,
    chains=2,
    prior=prior_list
))

## MAP estimate (this is the one Jo probably wants to use)
## "shortcut" Bayesian estimate: find the mode of the posterior
##  distribution
system.time(ff2 <- fitode(
    model=model,
    data=df,
    start=parms,
    prior=prior_list
))

## penalization based on prior (ignoring dmu/deta)
##  ignoring the change of scale (we're fitting some parameters
##  on the log scale)
system.time(ff3 <- fitode(
    model=model,
    data=df,
    start=parms,
    prior=prior_list,
    prior.density = FALSE
))

plot(ff, level=0.95)
plot(ff2, level=0.95)
plot(ff3, level=0.95)

summary(ff)
## Rhats are big, effective samples are small: this would be a problem
## if we were really doing this ...

lattice::xyplot(ff@mcmc) ## ugh.

coda::gelman.diag(ff@mcmc)

plot(ff@mcmc[,5])  ## plot chains and density of parameter 5 (size1)?
##  not sure whether this on link scale or original scale?

## QUESTIONS:
##  is the prior density always on the original scale? (i.e., how
##    is it connected with the link function for a parameter?)
##  how easy/hard is it to extend?
##  can we re-use information from the observation model stuff?
