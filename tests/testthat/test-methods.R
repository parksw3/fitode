stopifnot(require("testthat"), require("fitode"))

context("basic tests")
test_that("SI model", {
    harbin <- structure(list(week = 2:18,
                             Deaths = c(2, 7, 2, 6, 12, 68, 91,
                                        126, 229, 250, 212, 161, 101, 108, 46, 40, 14)),
                        .Names = c("week", "Deaths"),
                        row.names = c(NA, -17L), class = "data.frame")

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
                    par=c("beta", "gamma", "N", "i0")
    )

    start <- structure(c(1.85291753121151, 1.08599844275679, 2352.05157584658,
                         0.000480103245854955, 3.58305410314183),
                       .Names = c("beta", "gamma", "N", "i0", "ll.k"))

    ff <- fitode(Deaths~I,
        start=start,
        model=SI_model, loglik=select_model("nbinom"),
        data=harbin,
        tcol="week",
        links = list(
            beta="log",
            gamma="log",
            N="log",
            i0="logit"
        )
    )

    logLik(ff)


})
