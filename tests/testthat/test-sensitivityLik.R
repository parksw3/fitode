stopifnot(require("testthat"), require(numDeriv), require("fitode"))

test_that("SI model", {
    SI_model <- odemodel(
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
        unname(logLik.sensitivity(parms, SI_model, data))[-1],
        tolerance=1e-5
    )
})
