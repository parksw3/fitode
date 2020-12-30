stopifnot(require("testthat"), require(numDeriv), require("fitode"))

## context("log-likelihood tests") ## now deprecated
test_that("normal", {
    distrib <- select_model("dnorm")

    expect_equal(
        Eval(distrib, observation=1, mean=1, par=1),
        dnorm(1, 1, 1, log=TRUE),
        tolerance=1e-5
    )

    expect_equal(
        Eval(distrib, observation=2, mean=1, par=1),
        dnorm(2, 1, 1, log=TRUE),
        tolerance=1e-5
    )

    expect_equal(
        Eval(distrib, observation=1, mean=2, par=3),
        dnorm(1, 2, 3, log=TRUE),
        tolerance=1e-5
    )
})

test_that("poisson", {
    distrib <- select_model("dpois")

    expect_equal(
        Eval(distrib, observation=1, mean=1),
        dpois(1, 1, log=TRUE),
        tolerance=1e-5
    )

    expect_equal(
        Eval(distrib, observation=2, mean=1),
        dpois(2, 1, log=TRUE),
        tolerance=1e-5
    )

    expect_equal(
        Eval(distrib, observation=1, mean=3),
        dpois(1, 3, log=TRUE),
        tolerance=1e-5
    )
})

test_that("nbinom", {
    distrib <- select_model("dnbinom")

    expect_equal(
        Eval(distrib, observation=1, mean=1, par=2),
        dnbinom(1, mu=1, size=2, log=TRUE),
        tolerance=1e-5
    )

    expect_equal(
        Eval(distrib, observation=10, mean=1, par=0.1),
        dnbinom(10, mu=1, size=0.1, log=TRUE),
        tolerance=1e-5
    )

    expect_equal(
        Eval(distrib, observation=1, mean=10, par=1000),
        dnbinom(1, mu=10, size=1000, log=TRUE),
        tolerance=1e-5
    )
})

test_that("nbinom1", {
    distrib <- select_model("dnbinom1")

    dnbinom1 <- function(x, mu, phi, log) dnbinom(x, mu=mu, size=mu/phi, log=log)

    expect_equal(
        Eval(distrib, observation=1, mean=1, par=2),
        dnbinom1(1, mu=1, phi=2, log=TRUE),
        tolerance=1e-5
    )

    expect_equal(
        Eval(distrib, observation=10, mean=1, par=0.1),
        dnbinom1(10, mu=1, phi=0.1, log=TRUE),
        tolerance=1e-5
    )

    expect_equal(
        Eval(distrib, observation=1, mean=10, par=1000),
        dnbinom1(1, mu=10, phi=1000, log=TRUE),
        tolerance=1e-5
    )
})

test_that("gamma", {
    distrib <- select_model("dgamma")

    expect_equal(
        Eval(distrib, observation=1, mean=1, par=2),
        dgamma(1, shape=2, rate=2, log=TRUE),
        tolerance=1e-5
    )

    expect_equal(
        Eval(distrib, observation=5, mean=2, par=10),
        dgamma(5, shape=10, rate=5, log=TRUE),
        tolerance=1e-5
    )

    expect_equal(
        Eval(distrib, observation=100, mean=1, par=100),
        dgamma(100, shape=100, rate=100, log=TRUE),
        tolerance=1e-5
    )
})
