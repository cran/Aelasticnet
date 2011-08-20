\name{enetbqr}
\alias{enetbqr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Bayesian Adaptive Elastic-Net for high dimensional sparse binary quantile regression models.
}
\description{This function   introduces the idea of Bayesian variable selection and estimation in binary quantile regression models employing a likelihood function that is based on the
asymmetric Laplace distribution.}
\usage{
enetbqr(formula, tau = 0.5, burnin = 5000, mcmc = 25000,  keep = 1)
}
%- maybe also 'usage' for other objects documented here.

\arguments{
  \item{formula}{Model formula.}
  \item{tau}{The quantile of interest. Must be between 0 and 1.   }
  \item{burnin}{Number of Gibbs sampler iterations before output is saved.}
  \item{mcmc}{Length of desired Gibbs samples output.}
  \item{keep}{Thinning for chain (e.g. 1 means save every iteration).}
  %\item{...}{ Additional arguments to QRc.AL.}
}



\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{Rahim Alhamzawi and Keming Yu(2011). Gene selection using Bayesian adaptive elastic net binary quantile regression. Department of Mathematical sciences, Brunel University, Uxbridge UBB 3PH, UK.}

\author{Rahim Alhamzawi and Dries F. Benoit}
\note{
%%  ~~further notes~~
}

\examples{
###  Simulate data from the binary regression model.
n  = 200
k = 5
x = matrix(runif(n=n*k, min=0, max=10),n , k)
ystar = 1.5*x[, 1] + rnorm(n=n, mean=0, sd=1)
y <- as.numeric(ystar>0)

### Estimate the model
fit = enetbqr(y ~ 0+x, tau = 0.5)

### Check results
summary(fit$beta)
}
