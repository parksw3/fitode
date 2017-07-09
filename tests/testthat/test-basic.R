stopifnot(require("testthat"), require(numDeriv), require("fitode"))

context("basic tests")
test_that("SI model", {
    SI_model <- new("model.ode",
            name = "SI",
            model = list(
            S ~ - beta*S*I/N,
            I ~ beta*S*I/N - gamma*I
        ),
        state = c("S", "I"),
        par= c("beta", "gamma", "N")
    )
    ss <- solve(c(S=99, I=1), times = 1:10, SI_model, parms=c(beta=2, gamma=1, N=100))

    ff <- function(parms, y, model) {
        solve(y, times=1:10, model, parms=parms)@solution$I
    }

    expect_equal(
        numDeriv::jacobian(ff,c(beta=2, gamma=1, N=100), y=c(S=99, I=1), model=SI_model),
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
                    state = c("S", "E", "I"),
                    par= c("beta", "sigma", "gamma", "N")
    )

    ss <- solve(c(S=99, E=0, I=1), times = 1:10, SEI_model, parms=c(beta=2, sigma=1, gamma=1, N=100))

    expect_equal(
        numDeriv::jacobian(ff,c(beta=2, sigma=1, gamma=1, N=100), y=c(S=99, E=0, I=1), model=SEI_model),
        unname(as.matrix(ss@sensitivity$I)),
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
                     state = c("y1", "y2"),
                     par= c("mu")
    )

    ss <- solve(c(y1=2, y2=2), times = 1:10, VdP, parms=c(mu=0.1))

    ff <- function(parms, y, model) {
        solve(y, times=1:10, model, parms=parms)@solution$y1
    }

    expect_equal(
        numDeriv::jacobian(ff,c(mu=0.1), y=c(y1=2, y2=2), model=VdP)[,1],
        unlist(unname(ss@sensitivity$y1)),
        tolerance=1e-4
    )
})
