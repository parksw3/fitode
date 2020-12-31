stopifnot(require("testthat"), require(numDeriv), require("fitode"))

test_that("SI model", {
    SI_model <- odemodel(,
        name = "SI",
        model = list(
            S ~ - beta*S*I/N,
            I ~ beta*S*I/N - gamma*I
        ),
        observation = list(
            prev~dnbinom(mu=I, size=size)
        ),
        initial = list(
            S ~ N * (1 - i0),
            I ~ N * i0
        ),
        par= c("beta", "gamma", "N", "i0", "size")
    )

    ss <- ode.solve(SI_model, times = 1:10, parms=c(beta=2, gamma=1, N=1000, i0=0.001))

    ff <- function(parms) ode.solve(SI_model, times=1:10, parms=parms)@solution$I

    expect_equal(
        numDeriv::jacobian(ff, c(beta=2, gamma=1, N=1000, i0=0.001)),
        unname(as.matrix(ss@sensitivity$I[,1:4])),
        tolerance=1e-5
    )
})

test_that("SEI model", {
    SEI_model <- odemodel(
        name = "SEI",
        model = list(
            S ~ - beta*S*I/N,
            E ~ beta*S*I/N - sigma*E,
            I ~ sigma*E - gamma*I
        ),
        observation = list(
            prev~dnbinom(mu=I, size=size)
        ),
        initial = list(
            S ~ N * (1 - i0),
            E ~ 0,
            I ~ N * i0
        ),
        par= c("beta", "sigma", "gamma", "N", "i0", "size")
    )

    ss <- ode.solve(SEI_model, times = 1:10, parms=c(beta=2, sigma=1, gamma=1, N=1000, i0=0.001))

    ff <- function(parms) ode.solve(SEI_model, times=1:10, parms=parms)@solution$I

    expect_equal(
        numDeriv::jacobian(ff,c(beta=2, sigma=1, gamma=1, N=1000, i0=0.001)),
        unname(ss@sensitivity$I[,1:5]),
        tolerance=1e-5
    )

})

test_that("Van der Pol oscillator", {
    VdP <- odemodel("model.ode",
        name = "VdP",
        model = list(
            y1 ~ y2,
            y2 ~ mu*(1-y1^2)*y2-y1
        ),
        observation = list(
            obs1 ~ dnorm(mean=y1, sd=sigma)
        ),
        initial = list(
            y1 ~ 2,
            y2 ~ 2
        ),
        par= c("mu", "sigma")
    )

    ss <- ode.solve(VdP, times = 1:10, parms=c(mu=0.1))

    ff <- function(parms) ode.solve(VdP, times=1:10, parms=parms)@solution$y2

    expect_equal(
        c(numDeriv::jacobian(ff, c(mu=0.1))),
        unname(ss@sensitivity$y2[,1]),
        tolerance=1e-3
    )
})
