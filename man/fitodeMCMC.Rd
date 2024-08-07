% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fitodeMCMC.R
\name{fitodeMCMC}
\alias{fitodeMCMC}
\title{Fit ordinary differential equations model using MCMC}
\usage{
fitodeMCMC(
  model,
  data,
  start,
  tcol = "times",
  proposal.vcov,
  prior = list(),
  chains = 1,
  iter = 2000,
  burnin = iter/2,
  thin = 1,
  refresh = max(iter/10, 1),
  prior.only = FALSE,
  link,
  fixed = list(),
  solver.opts = list(method = "rk4"),
  solver = ode,
  ...
)
}
\arguments{
\item{model}{ode model}

\item{data}{data frame with time column and observation column}

\item{start}{named vector of starting parameter values}

\item{tcol}{time column}

\item{proposal.vcov}{variance-covariance matrix of a multivariate normal proposal distribution}

\item{prior}{list of formulas specifying prior distributions}

\item{chains}{(numeric) number of chains}

\item{iter}{(numeric) number of iterations per chain}

\item{burnin}{(numeric) number of burnin interations}

\item{thin}{(numeric) thining interval between consecutive observations}

\item{refresh}{(numeric) refresh interval}

\item{prior.only}{(logical) sample from prior distribution only?}

\item{link}{named vector or list of link functions for model parameters}

\item{fixed}{named vector or list of model parameters to fix and their values}

\item{solver.opts}{options for ode integration. See \code{\link{ode}}}

\item{solver}{ode solver}

\item{...}{additional arguments (unused)}
}
\value{
An object of class ``fitodeMCMC'' as described in \code{\link{fitodeMCMC-class}}.
}
\description{
This function fits ordinary differential equations models to a uni- or
multi-variate time series by MCMC using the Metropolis-Hastings update
rule. It searches through the parameter space
on link scales, which can provide more efficient posterior sampling.
}
\seealso{
\code{\link{fitodeMCMC}}
}
