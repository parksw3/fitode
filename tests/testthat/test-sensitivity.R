stopifnot(require("testthat"), require(numDeriv), require("fitode"))

context("basic tests")
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
    ss <- solve(SI_model, times = 1:10, parms=c(beta=2, gamma=1, N=1000, i0=0.001))

    ff <- function(parms, model, var) {
        solve(model, times=1:10, parms=parms)@solution[[var]]
    }

    expect_equal(
        numDeriv::jacobian(ff, c(beta=2, gamma=1, N=1000, i0=0.001), model=SI_model, var="I"),
        unname(as.matrix(ss@sensitivity$I)),
        tolerance=1e-5
    )
})

test_that("SEI model", {
    SEI_model <- new("model.ode",
        name = "SEI",
        model = list(
            S ~ - beta*S*I/N,
            E ~ beta*S*I/N - sigma*E,
            I ~ sigma*E - gamma*I
        ),
        initial = list(
            S ~ N * (1 - i0),
            E ~ 0,
            I ~ N * i0
        ),
        par= c("beta", "sigma", "gamma", "N", "i0")
    )

    ss <- solve(SEI_model, times = 1:10, parms=c(beta=2, sigma=1, gamma=1, N=1000, i0=0.001))

    expect_equal(
        numDeriv::jacobian(ff,c(beta=2, sigma=1, gamma=1, N=1000, i0=0.001), model=SEI_model, var="I"),
        unname(ss@sensitivity$I),
        tolerance=1e-5
    )

})

test_that("Van der Pol oscillator", {
    VdP <- new("model.ode",
        name = "VdP",
        model = list(
            y1 ~ y2,
            y2 ~ mu*(1-y1^2)*y2-y1
        ),
        initial = list(
            y1 ~ 2,
            y2 ~ 2
        ),
        par= "mu"
    )

    ss <- solve(VdP, times = 1:10, parms=c(mu=0.1))

    expect_equal(
        numDeriv::jacobian(ff, c(mu=0.1), model=VdP, var="y2"),
        unname(ss@sensitivity$y2),
        tolerance=1e-3
    )
})
