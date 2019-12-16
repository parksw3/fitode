library(fitode)

## making fake data
set.seed(123)
true.m <- 0.5
true.A0 <- 200
true.sd <- 15

exp.data <- data.frame(
    x=rep(0:4, 10),
    y=rnorm(5*10, true.A0 * exp(-true.m * rep(0:4, 3)), sd=true.sd)
)

plot(exp.data)

## how we used to set up gaussian models
exp.model <- new("model.ode",
                 name = "exp",
                 model = list(
                     A ~ -m * A
                 ),
                 observation = list(
                     y ~ dnorm(mean=A, sd=sd)
                 ),
                 initial = list(
                     A ~ A0
                 ),
                 par=c("m", "A0", "sd")
)

## now we don't have to estimate sd (using dnorm2 instead of dnorm)
exp.model2 <- new("model.ode",
                 name = "exp",
                 model = list(
                     A ~ -m * A
                 ),
                 observation = list(
                     y ~ dnorm2(mean=A)
                 ),
                 initial = list(
                     A ~ A0
                 ),
                 par=c("m", "A0")
)

## comparison of fits:
exp.fit <- fitode(
    exp.model,
    exp.data,
    start=c(m=0.5, A0=200, sd=15),
    tcol="x"
)

exp.fit2 <- fitode(
    exp.model2,
    exp.data,
    start=c(m=0.5, A0=200),
    tcol="x"
)
## slightly different
coef(exp.fit)
coef(exp.fit2)

## estimated sd
coef(exp.fit)[3]
## profiled sd
sqrt(sum((exp.data$y-rep(predict(exp.fit2)$y[,2], 10))^2)/49)
## dnorm2 actually uses the unbiased estimator (dividing by n-1) for sd rather than MLE (dividing by n)
## that's why they're slightly different
