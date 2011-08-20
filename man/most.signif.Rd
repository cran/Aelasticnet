\name{most.signif}
\alias{most.signif}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Selection the significant variables}
\description{\code{most.signif} This function identify the most significant  variables that contribute most to the classification.}
\usage{
most.signif(sv, tau = 0.5, burnin = 3000, mcmc = 10000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sv}{Number of significant variables.}
  \item{tau}{The quantile of interest. Must be between 0 and 1.   }
  \item{burnin}{Number of Gibbs sampler iterations before output is saved.}
  \item{mcmc}{Length of desired Gibbs samples output.}
  %\item{keep}{Thinning for chain (1 = save every iteration).}
  %\item{...}{ Additional arguments to QRc.AL.}
}
\examples{
###  Simulation data.
## Example: Most significant variables.
n  = 200
k = 5
x = matrix(runif(n=n*k, min=0, max=10),n , k)
ystar = 1.5*x[, 5] + rnorm(n=n, mean=0, sd=1)
y <- as.numeric(ystar>0)

most.signif(sv = 3, tau = 0.5)
}