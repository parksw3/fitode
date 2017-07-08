stopifnot(require("testthat"), require("fitode"))

context("basic tests")
test_that("SI model", {
    library(deSolve) ## TODO: fix this... put it in the imports
    library(numDeriv) ## get rid of this some time later...
    SI_model <- new("de-model",
            name = "SI",
            model = list(
            S ~ - beta*S*I/N,
            I ~ beta*S*I/N - gamma*I
        ),
        state = c("S", "I"),
        par= c("beta", "gamma", "N")
    )
    ss <- solve(c(S=99, I=1), times = 1:10, SI_model, parms=c(beta=2, gamma=1, N=100))

    ff <- function(parms) {
        solve(c(S=99, I=1), times=1:10, SI_model, parms=parms)[,2]
    }

    expect_equal(
        numDeriv::jacobian(ff,c(beta=2, gamma=1, N=100)),
        unname(ss[,4:6]),
        tolerance=1e-5
    )
})
