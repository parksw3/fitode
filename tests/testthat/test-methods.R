stopifnot(require("testthat"), require("fitode"))

test_that("SI model", {
    harbin <- structure(list(week = 2:18,
        Deaths = c(2, 7, 2, 6, 12, 68, 91,
            126, 229, 250, 212, 161, 101, 108, 46, 40, 14)),
        .Names = c("week", "Deaths"),
        row.names = c(NA, -17L),
        class = "data.frame")

    SI_model <- odemodel(
        name = "SI",
        model = list(
            S ~ - beta*S*I/N,
            I ~ beta*S*I/N - gamma*I
        ),
        observation = list(
            Deaths ~ dnbinom(mu=gamma*I, size=size)
        ),
        initial = list(
            S ~ N * (1 - i0),
            I ~ N * i0
        ),
        par=c("beta", "gamma", "N", "i0", "size"),
        link=c(i0="logit")
    )

    start <- c(beta=2, gamma=1, N=1e4, i0=1e-3, size=10)

    suppressWarnings(ff <- fitode(
        SI_model,
        start=start,
        data=harbin,
        tcol="week",
        quietly = TRUE
        
    ))

    expect_equal(
        coef(ff),
        structure(c(1.85082710742599, 1.08315701471954, 2162.03690454815,
                    0.000482746970936824, 35.9813684389218),
                  .Names = c("beta", "gamma", "N", "i0", "size")),
        tolerance = 1e-6)

    expect_equal(
        coef(ff, "links"),
        structure(c(0.615632624272599, 0.0798799387838958, 7.67880606768721,
            -7.63553504774462, 3.58300126112103),
            .Names = c("log.beta", "log.gamma", "log.N", "logit.i0", "log.size")),
        tolerance = 1e-6)

    expect_equal(
        logLik(ff),
        -67.6503511573217,
        tolerance = 1e-6
    )


    expect_equal(
        tolerance = 1e-6,
        predict(ff, level=0.95)[[1]],
        structure(
            list(times = 2:18,
                mean = c(1.13050913725776, 2.42831844773126,
                    5.19681134127028, 11.0346210723243, 23.0460972694204, 46.5271091921174,
                    87.9263920319632, 147.842264944473, 208.456795463428, 237.435011324319,
                    220.45598366043, 174.414692008427, 123.684857749301, 81.8307738696627,
                    51.8905997500834, 32.0750023485528, 19.5258283338705),
                `2.5 %` = c(0.376399546733032,1.10052485552639, 2.95578804185289, 7.4424083345586, 17.4742055021702,
                    37.232745585215, 70.3978733348157, 118.255855900933, 171.28077590353,
                    198.159013874503, 180.609882060902, 141.530263861177, 101.890963842818,
                    68.15590115031, 41.7836369098972, 23.633131796506, 12.5932104646235),
                `97.5 %` = c(1.88461872778249, 3.75611203993613, 7.43783464068767,
                    14.6268338100899, 28.6179890366706, 55.8214727990198, 105.454910729111,
                    177.428673988013, 245.632815023326, 276.711008774134, 260.302085259957,
                    207.299120155678, 145.478751655784, 95.5056465890154, 61.9975625902696,
                    40.5168729005996, 26.4584462031175)),
            .Names = c("times", "estimate", "2.5 %", "97.5 %"), row.names = c(NA, -17L),
            class = "data.frame")
    )
})
