library(bbmle)
library(emdbook)
library(MCMCpack)
library(ggplot2); theme_set(theme_bw())
library(gridExtra)
set.seed(101)
x <- rnorm(25)

m <- mle2(x~dnorm(mean=mean, sd=exp(ls)),
     start=list(mean=0, ls=0),
     data=data.frame(x=x))

loglik <- function(x, y, data) {
    -sum(dnorm(data, x, exp(y), log=TRUE))
}

xlim <- seq(-0.6, 0.5, by=0.02)
ylim <- seq(-0.5, 0.5, by=0.02)

loglik_surface <- apply2d(loglik,
        x=xlim,
        y=ylim,
        data=x)

pp <- profile(m)

m2 <- MCMCmetrop1R(function(param, data) sum(dnorm(data, param[1], exp(param[2]), log=TRUE)),
             theta.init=c(0, 0),
             burnin=50000,
             mcmc=200000,
             thin=10,
             data=x,
             V=vcov(m))

surfacedf <- data.frame(
    x=rep(xlim, length(ylim)),
    y=rep(ylim, each=length(xlim)),
    z=round(c(loglik_surface))
)

profiledf <- as.data.frame(pp@profile$mean[,2])

mledf <- data.frame(
    mean=coef(m)[1],
    ls=coef(m)[2]
)

margquant <- data.frame(
    lwr=quantile(m2[,1], 0.025),
    upr=quantile(m2[,1], 0.975)
)

prof.conf <- confint(pp)[1,]

t.conf <- data.frame(
    conf=t.test(x)$conf.int,
    y=-0.4
)

z.conf <- data.frame(
    conf=mean(x) + c(-1, 1) * 1.96 * sqrt(vcov(m)[1]),
    y=-0.45
)

gplot <- ggplot(surfacedf) +
    geom_raster(aes(x, y, fill=z)) +
    geom_line(data=profiledf, aes(mean, ls), lty="longdash", lwd=1) +
    geom_point(data=mledf, aes(mean, ls), size=10) +
    geom_vline(xintercept=margquant$lwr, lty="dotted", col=2) +
    geom_vline(xintercept=margquant$upr, lty="dotted", col=2) +
    geom_vline(xintercept=prof.conf[1], lty="longdash") +
    geom_vline(xintercept=prof.conf[2], lty="longdash") +
    geom_line(data=t.conf, aes(conf, y)) +
    annotate("text", x=coef(m)[1], y=-0.37, label="t-test") +
    geom_line(data=z.conf, aes(conf, y), col="blue") +
    annotate("text", x=coef(m)[1], y=-0.48, label="z-score (with estimated s.e.)", col="blue") +
    scale_x_continuous(expression(mu), expand=c(0,-0.1), limits=c(-0.65, 0.55)) +
    scale_y_continuous(expression(log(sigma)), expand=c(0,0)) +
    scale_fill_gradient(low="white", high="black") +
    theme(
        legend.position = "none"
    )

ghist <- ggplot(data.frame(x=m2[,1])) +
    geom_histogram(aes(var1), bins=30, fill=NA, col="black") +
    geom_vline(xintercept=margquant$lwr, lty="dotted", col=2) +
    geom_vline(xintercept=margquant$upr, lty="dotted", col=2) +
    scale_x_continuous(expand=c(0,-0.1), limits=c(-0.65, 0.55)) +
    scale_y_continuous("") +
    theme(
        panel.background = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_line(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(margin = margin(t = 0, r = 37, b = 0, l = 0)),
        axis.ticks.y = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank()
    )

gfinal <- arrangeGrob(ghist, gplot,
    nrow=2, heights=c(0.3, 0.7)
)

ggsave("normal.pdf", gfinal, width=8, height=6)
