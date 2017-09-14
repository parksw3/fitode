stopifnot(require("testthat"), require(numDeriv), require("fitode"))

context("sinusoidal SI model tests")
test_that("sinusoidal model", {
    ## based on sinusoidal model in EpiDynamics

    SI_nm_model <- new("model.ode",
        name = "SI",
        model = list(
            S ~ mu-beta*S*I - mu * S,
            I ~ beta*S*I-mu*I-gamma*I
        ),
        initial = list(
            S ~ S0,
            I ~ I0
        ),
        par=c("beta", "mu", "gamma", "S0", "I0")
    )

    sinusoidal_model <- Transform(
        SI_nm_model,
        list(beta~beta0 * (1 + beta1 * sin(omega * t))),
        par=c("beta0", "beta1" ,"omega", "gamma", "mu", "S0", "I0")
    )

    t <- 0:(10*365)

    parms <- c(beta0=17/13, beta1=0.1, omega=2*pi/365, gamma=1/13, mu=1/(50*365), S0=1/17, I0=1e-4)

    ss <- ode.solve(sinusoidal_model, t, parms)

    ff <- function(parms) ode.solve(sinusoidal_model,
                                    times=t,
                                    parms=parms,
                                    keep_sensitivity=FALSE)@solution$I

    expect_equal(
        numDeriv::jacobian(ff, parms),
        unname(as.matrix(ss@sensitivity$I)),
        tol=1e-4
    )

})
