\name{enetqr}
\alias{enetqr}
\alias{Model.selction}
\title{Bayesian Adaptive Elastic-Net for high dimensional sparse quantile regression models}

\description{This function   introduces the idea of Bayesian variable selection and estimation in  quantile regression models employing a likelihood function that is based on the
asymmetric Laplace distribution. A novel aspect of the Bayesian
adaptive elastic net quantile regression is to treat all  the hyperparameters as unknowns and let the data estimate them along with other parameters. We also treat the normalizing constant as unknown.
To improve the speed of the routine, the Markov Chain Monte Carlo (MCMC) part of the algorithm is programmed in FORTRAN and is called from within the R function enetqr.
}
\usage{
enetqr(formula, tau=0.5, burnin = 10000, mcmc=25000,keep=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{Model formula.}
  \item{tau}{The quantile of interest. Must be between 0 and 1.   }
  \item{burnin}{Number of Gibbs sampler iterations before output is saved.}
  \item{mcmc}{Length of desired Gibbs samples output.}
  \item{keep}{Thinning for chain (E.g. 1 means save every iteration).}
  %\item{...}{ Additional arguments to QRc.AL.}
}
\details{}
%%\value{A list containing:
%%  ~Describe the value returned
%%  If it is a LIST, use
 %% \item{betadraw} { mcmc x nvar(X) array of beta draws}
 %% \item{sdraw }{mcmc x nvar(X) array of beta draws}
%% ...
%%}
\references{Rahim Alhamzawi and Keming Yu(2011). Gene selection using Bayesian adaptive elastic net binary quantile regression. Department of Mathematical sciences, Brunel University, Uxbridge UBB 3PH, UK.}

\author{Rahim Alhamzawi and   Dries F. Benoit }

\note{Please contact Rahim Alhamzawi <rahim.al_hamzawi@brunel.ac.uk> with comments or suggestions.}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ Q. Li,  R. Xi,  and N. Lin.  Bayesian regularized quantile regression.  \emph{Bayesian Analysis},  5(3):533--56,  2010.

 Roger Koenker (2011). quantreg: Quantile Regression. R package version 4.71. http://CRAN.R-project.org/package=quantreg.

The Base package in R.

}
\examples{
###  Real data from heteroskedastic regression..
## Bayesian estimation
library(quantreg)
data(engel)
attach(engel)
y=foodexp
x=income
# comparing between rq and enetqr.
# tau=0.05
fit = enetqr(foodexp ~ income, tau = 0.05)
summary(fit$beta)
# credible interval
round(apply(fit$beta,2,quantile,c(0.025,0.5,0.975)),3)
summary(rq(foodexp ~ income, tau = 0.05))

# tau=0.50
fit = enetqr(foodexp ~ income, tau = 0.50)
summary(fit$beta)
# credible interval
round(apply(fit$beta,2,quantile,c(0.025,0.5,0.975)),3)
summary(rq(foodexp ~ income, tau = 0.50))

# tau=0.95
fit = enetqr(foodexp ~ income, tau = 0.95)
summary(fit$beta)
# credible interval
round(apply(fit$beta,2,quantile,c(0.025,0.5,0.975)),3)
summary(rq(foodexp ~ income, tau = 0.95))

# Initiate plot
## Plot datapoints
plot(income, foodexp, main="", cex=.6, xlab="income")

# Write loop to analyze 5 quantiles
for (i in 1:5) {
if (i==1) p = .05
if (i==2) p = .25
if (i==3) p = .50
if (i==4) p = .75
if (i==5) p = .95

fit = enetqr(foodexp~income,tau=p)
#print(summary(fit$beta))

abline(a=mean(fit$beta[,1]),b=mean(fit$beta[,2]),lty=i,col=i)}
# Estimate and plot OLS model
outOLS = lm(foodexp~income)
abline(outOLS,lty=1,lwd=2,col=6)

# Add legend to plot
legend(x=0,y=max(y),legend=c(.05,.25,.50,.75,.95,"OLS"),lty=c(1,2,3,4,5,1),
lwd=c(1,1,1,1,1,2),col=c(1:6),title="Quantile")


}


