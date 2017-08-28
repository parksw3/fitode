stopifnot(require("testthat"), require(numDeriv), require("fitode"))

context("likelihood sensitivity tests")
test_that("SI model", {
    SI_model <- new("model.ode",
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

    parms <- c(beta=2, gamma=1, N=1000, i0=0.001)
    Deaths <- c(1:10)
    times <- c(1:10)

    ff <- function(parms, model) {
        ss <- ode.solve(model, times, parms=parms)@solution
        -sum(dpois(Deaths, parms[1] * ss$S * ss$I, log=TRUE))
    }

    expect_equal(
        as.vector(numDeriv::jacobian(ff, parms, model=SI_model)),
        unname(fitode:::logLik.sensitivity(parms, Deaths~beta*S*I, SI_model, select_model("poisson"), Deaths, times)[-1]),
        tolerance=1e-5
    )
})
