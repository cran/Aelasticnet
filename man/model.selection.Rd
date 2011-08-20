\name{model.selection}
\alias{model.selection}
\title{Bayesian adaptive Elastic-Net for Sparse Estimation}
\description{\code{model.selection} This function uses  Bayesian adaptive Elastic-Net to choose the optimal model for continuous data. }
\usage{
model.selection(formula, tau = 0.5, burnin = 10000, mcmc = 25000, keep = 1)
}
\arguments{
  \item{formula}{Model formula.}
  \item{tau}{The quantile of interest. Must be between 0 and 1.   }
  \item{burnin}{Number of Gibbs sampler iterations before output is saved.}
  \item{mcmc}{Length of desired Gibbs samples output.}
  \item{keep}{Thinning for chain (1 = save every iteration).}
  %\item{...}{ Additional arguments to QRc.AL.}
}
\examples{
## Example: Model selection for continuous data.

n=200
k=15
x=matrix(runif(n*k,0,10),n,k)
x=cbind(1,x)
y=10-0.5*x[,2]+1/11*(11+x[,2])*rnorm(n,0,1)

p=0.75
model.selection(y~x,tau=p, burnin = 1000, mcmc = 10000)
library(quantreg)
rq(y~0+x,tau=p)$coeff

p=0.5
model.selection(y~x,tau=p, burnin = 1000, mcmc = 10000)
rq(y~0+x,tau=p)$coeff

p=0.25
model.selection(y~x,tau=p, burnin = 1000, mcmc = 10000)
rq(y~0+x,tau=p)$coeff
}





