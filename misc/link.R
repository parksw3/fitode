library(fitode)

model <- new("model.ode",
    name = "SI",
    model = list(
        S ~ - beta*S*I/N,
        I ~ beta*S*I/N - gamma*I
    ),
    initial = list(
        S ~ N * (1 - i0),
        I ~ N * i0
    ),
    par= c("beta", "gamma", "N", "i0")
)

linkfun <- function(link=c("log", "logit")) {
    link <- match.arg(link)
    switch(link,
        log={
            list(
                transform=substitute(exp(x)),
                inverse=substitute(log(x))
            )
        },
        logit={
            list(
                transform=substitute(1/(1+exp(-x))),
                inverse=substitute(log(x)-log(1-x))
            )
        }
    )
}

links <- list(
    beta="log",
    N="log"
)

formula <- Deaths~beta*S*I/N

start <- c(beta=2, gamma=1, N=1e5, i0=1e-4, ll.k=2)

transform <- vector('list', length(links))
inverse <- vector('list', length(links))

newpar <- oldpar <- model@par

for(i in 1:length(links)) {
    par <- names(links)[i]
    link <- links[[i]]

    tpar <- paste(link, par, sep=".")

    newpar[which(newpar==par)] <- tpar

    m <- Map(subst, e=linkfun(link), transforms=list(list(x=as.name(tpar)), list(x=as.name(par))))

    transform[[i]] <- to.formula(par, m$transform)
    inverse[[i]] <- to.formula(tpar, m$inverse)
}

model <- Transform(model, transforms=transform, par=newpar)

newformula <- subst(formula[[3]], trans(transform, oldpar))
formula <- to.formula(formula[[2]], newformula)

odestart <- start[oldpar]

## TODO: re-write this section ...

ss <- lapply(inverse, subst, as.list(odestart))
ss <- lapply(ss, as.formula)
tss <- trans(ss, newpar)
newstart <- as.list(odestart)
names(newstart) <- newpar
newstart[names(tss)] <- tss
start <- c(sapply(newstart, eval), start[loglik@par])


## here:

start
transform
inverse

transpar <- function(parms, newname, inverse) {
    parlist <- lapply(newname, function(x) as.name(x))
    newpar <- lapply(lapply(parlist, subst, trans(inverse,newname)), subst, as.list(parms))
    newpar <- sapply(newpar, eval)
    names(newpar) <- newname
    newpar
}

transpar(transpar(start, newpar, inverse), oldpar, transform)

parms <- start
inverse <- inverse

parms
inverse

