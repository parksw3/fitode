stopifnot(require("testthat"), require(numDeriv), require("fitode"))

context("likelihood sensitivity tests")
test_that("SI model", {
    SI_model <- new("model.ode",
        name = "SI",
        model = list(
            S ~ - beta*S*I/N,
            I ~ beta*S*I/N - gamma*I,
            C ~ gamma * I
        ),
        initial = list(
            S ~ N * (1 - i0),
            I ~ N * i0,
            C ~ 0
        ),
        observation = list(
            Deaths ~ dpois(lambda=C)
        ),
        diffnames="C",
        par= c("beta", "gamma", "N", "i0")
    )

    parms <- c(beta=2, gamma=1, N=1000, i0=0.001)
    Deaths <- c(NA, 1:10)
    times <- c(0:10)
    data <- data.frame(times=times, Deaths=Deaths)

    ff <- function(parms, model) {
        ss <- ode.solve(model, times, parms=parms)@solution
        -sum(dpois(Deaths, ss$C, log=TRUE), na.rm=TRUE)
    }

    expect_equal(
        as.vector(numDeriv::jacobian(ff, parms, model=SI_model)),
        unname(fitode:::logLik.sensitivity(parms, SI_model, data, fixed=NULL))[-1],
        tolerance=1e-5
    )
})

test_that("dnorm2", {
    model <- new("model.ode",
                    name = "constant",
                    model = list(
                        X ~ k
                    ),
                    initial = list(
                        X ~ 0
                    ),
                    observation = list(
                        Y ~ dnorm2(mean=X)
                    ),
                    par= c("k")
    )

    parms <- c(k=2)
    data <- data.frame(times=1:3, Y=c(0, 1, 5))

    ff <- function(parms, model) {
        ss <- ode.solve(model, times, parms=parms)@solution
        -sum(dpois(Deaths, ss$C, log=TRUE), na.rm=TRUE)
    }

    expect_equal(
        as.vector(numDeriv::jacobian(ff, parms, model=SI_model)),
        unname(fitode:::logLik.sensitivity(parms, model, data, fixed=NULL))[-1],
        tolerance=1e-5
    )
})
