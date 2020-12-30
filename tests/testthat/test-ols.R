stopifnot(require("testthat"), require(numDeriv), require("fitode"))

test_that("compare ols and dnorm", {
    ## UGH, I hate testthat() scoping model.  Putting this all inside a single test_that() clause
    ## to avoid headaches
    lin_model <- odemodel(
        name="linear trend",
        model=list(
            x ~ ax,
            y ~ ay
        ),
        observation=list(
            x_obs ~ ols(mean=x),
            y_obs ~ ols(mean=y)
        ),
        initial=list(
            x ~ 1,
            y ~ 1
        ),
        link=c(ax="identity", ay="identity"),
        par=c("ax","ay")
    )
    set.seed(101)
    nt <- 20
    dd <- data.frame(times=1:nt)
    dd$x_obs <- rnorm(nt, dd$times, sd=1)
    dd$y_obs <- rnorm(nt, 2*dd$times, sd=1)
    ## works fine in global environment
    f1 <- fitode(lin_model,
                 data=dd,
                 start=c(ax=1,ay=2),
                 quietly = TRUE)
    f2 <- update(f1,
                 observation=list(
                     x_obs ~ dnorm(mean=x,sd=sd1),
                     y_obs ~ dnorm(mean=y,sd=sd1)
                 ),
                 par=c("ax","ay","sd1"),
                 start=c(ax=1,ay=2,sd1=1)
                 )
    ## test_that("fit equals reference values",
    expect_equal(coef(f1), c(ax = 0.974617382010125, ay = 2.07757762305141))
    ## test_that("ols fit equals equal-var dnorm fit",
    expect_equal(coef(f1), coef(f2)[1:2])
    ## FIXME: should be closer?
    ## sqrt(2*nt/(2*nt-1)) should account for var() computing RSS/(n-1) rather than RSS/n
    ## test_that("correction to vcov()", {
    expect_equal(stdEr(f1)*sqrt(2*nt/(2*nt-1)), stdEr(f2)[-3], tolerance=1e-2)
})

## FIXME: construct an lm() that fits the same model?
