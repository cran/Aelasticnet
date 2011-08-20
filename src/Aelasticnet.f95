! Written by 
! Rahim Al-Hamzawi
! Department of mathematics
! Brunel University - UK
! and
! Dries F. Benoit 
! Faculty of economics and business administration
! Ghent University - BELGIUM

! MCMC sampler for Bayesian adaptive elastic net quantile regression. 
! This code is based on the following paper:
! Rahim Alhamzawi and Keming Yu(2011). Gene selection using Bayesian adaptive elastic net binary quantile regression.
! Department of Mathematical sciences, Brunel University, Uxbridge UBB 3PH, UK.

! Input arguments:
!	- n		: number of units of analysis
!	- k		: number of independent variables
!	- r		: number of MCMC iterations
!	- keep		: thinning parameter of MCMC draws
!	- y		: dependent variable
!	- p		: quantile of interest
!	- step_delta	: Metropolis-Hastings stepsize for delta
!	- x		: matrix of regressors (1:n, 1:k)
!       - sig_shape     : shape hyperparameter for sigma
!       - sig_rate      : rate hyperparameter for sigma

! Output arguments:
!	- betadraw	: the matrix of regression parameter estimates
!       - lambda12draw   : the matrix of lambda12draws
!       - deltadraw     : the vector of delta's
!       - taudraw       : the vector of tau's
!       - rejrate       : the rejection rate of the M-H step

SUBROUTINE Aelasticnet (n, k, r, keep, y, p, step_delta, x, sig_shape, sig_rate, &
                        betadraw, lambda12draw, deltadraw, taudraw, rejrate)

IMPLICIT NONE

! Precision statement:
INTEGER, PARAMETER :: dp = KIND(1.0d0)

! Input arguments:
INTEGER, INTENT(IN) :: n, k, r, keep
REAL(dp), INTENT(IN) :: p, sig_shape, sig_rate, step_delta
REAL(dp), INTENT(IN), DIMENSION(1:n) :: y
REAL(dp), INTENT(IN), DIMENSION(1:n,1:k) :: x

! Output arguments:
REAL(dp), INTENT(OUT) :: rejrate
REAL(dp), INTENT(OUT), DIMENSION(1:(r/keep)) :: deltadraw, taudraw
REAL(dp), INTENT(OUT), DIMENSION(1:(r/keep),1:k) :: betadraw, lambda12draw

! Internal arguments:
INTEGER :: naccept, i1, i2, i3, cnt
INTEGER, DIMENSION(1:(k-1)) :: selec
REAL(dp) :: sigma, delta, tau, theta, phisq,  tlambda1, &
            tvar, tmean, rnor, tshape, trate1, deltanew, llnew, llold, lldiff, lambda2new,nu2, nu12,&
            alpha, unif,lambda2
REAL(dp), DIMENSION(1:n) :: z, tmu1
REAL(dp), DIMENSION(1:k) :: beta, s, lambda12, trate2, tmu2, tlambda2

LOGICAL :: ISNAN


! -- SET STARTING VALUES
z = 1.0_dp
beta = 1.0_dp
s = 1.0_dp
lambda12 = 1.0_dp
lambda2=1.0_dp
nu2=1.0_dp
nu12=1.0_dp
sigma = 1.0_dp
delta = 1.0_dp
tau = 1.0_dp
naccept = 0


! -- CALCULATE USEFUL QUANTITIES
!if(p .GT. 0.95_dp) THEN p=0.96_dp
!if(p .LT. 0.05_dp) THEN p=0.04_dp
theta = 1.0_dp - 2.0_dp*p
phisq = p*(1.0_dp-p)


! -- START OF MCMC CHAIN
DO i1 = 1,r

  ! draw the latent z
  tmu1 = sqrt(1/(y-matmul(x, beta))**2.0_dp)
  tlambda1 = sigma/2
  DO i2 = 1,n
    CALL rInvGaussian(tmu1(i2), tlambda1, z(i2)) 
  END DO
  z = z**(-1.0_dp)
  

  ! draw mixing variable s
 
    tlambda2 = beta ** 2.0_dp + 1e-6
    tmu2 = sqrt(tlambda2/(sigma * lambda12))
    
  DO i2 = 1,k
    CALL rInvGaussian(tmu2(i2), tlambda2(i2), s(i2))
  END DO
s = s**(-1.0_dp)
  ! draw the regression parameters beta
  DO i2 = 1,k
  cnt = 0
    DO i3 = 1,k
      IF (i2 .NE. i3) THEN
        cnt = cnt + 1
        selec(cnt) = i3
      END IF
    END DO

    tvar = (sigma*sum(x(1:n,i2)**2.0_dp/(2*z)) + &
           s(i2)+nu12/nu2*lambda2)**(-1.0_dp)
    tmean = tvar*sigma*sum(x(1:n,i2)/(2*z) * &
            (y - matmul(x(1:n,selec),beta(selec)) - theta*z))
    CALL rnorm(rnor)
    beta(i2) = tmean + rnor*sqrt(tvar)
    
if (isnan(beta(i2))) beta(i2) =  betadraw((i1-1/keep),i2)
if (s(i2) .GT. 3.5e3_dp) beta(i2)=0.0_dp

  END DO

  ! draw sigma
  tshape = sig_shape  + 3.0_dp/2.0_dp*real(n,dp)
  trate1 = sum(((y - matmul(x,beta) - theta*z)**2.0_dp/(2.0_dp*z)) + &
         phisq*z)  + sig_rate
  CALL rgamma(tshape, trate1, sigma)

  sigma=sqrt(sigma)

   IF (sigma .LT. 1.0_dp) THEN
     sigma = sigma**(-1.0_dp)
   END IF

  ! draw lambda12
  tshape = 1.0_dp + delta
  trate2 = nu12/(2.0_dp*s*nu2)+ tau
  DO i2 = 1,k
    CALL rgamma(tshape, trate2(i2), lambda12(i2))
  END DO
 

! draw nu2
 tshape = sig_shape + real(k,dp) 
 trate1  = sum((nu12*lambda12)/(2*s) + nu12*lambda2*beta**2/2)  + sig_rate
  CALL rgamma(tshape, trate1**(-1.0_dp), nu2)
  nu2=nu2**(-1.0_dp)

! draw nu12
 tshape = sig_shape + real(k,dp) + 0.5_dp
 trate1 = sum(lambda12/(2*nu2*s) + lambda2*beta**2/(2*nu2))  + sig_rate
  CALL rgamma(tshape, trate1, nu12)

  ! draw the shrinkage parameter  LambdA2: Metropolis-Hastings
  LambdA2new = -1.0_dp
  DO WHILE (LambdA2new .LE. 0.0_dp)
    CALL rnorm(rnor)
    lambdA2new = lambda2 + rnor*step_delta
  END DO

  CALL lambd(lambdA2new,k, nu12, nu2, beta, llnew)
  CALL lambd(lambda2, k, nu12, nu2, beta, llold)

  lldiff = llnew - llold
  alpha = min(1.0_dp, exp(lldiff))

  IF (alpha .LT. 1.0_dp) THEN
    CALL random_number(unif)
  ELSE
    unif = 0.0_dp
  END IF

  IF (unif .LE. alpha) THEN
    lambda2 = LambdA2new
    naccept = naccept + 1
  END IF
  

  ! draw of hyperparameter tau
  tshape = delta*real(k,dp)
  trate1 = sum(lambda12)
  CALL rgamma(tshape, trate1, tau)
  if (tau .GT. 0.01) tau=0.01
   

  ! draw of shape hyperparameter delta: Metropolis-Hastings
  deltanew = -1.0_dp
  DO WHILE (deltanew .LE. 0.0_dp)
    CALL rnorm(rnor)
    deltanew = delta + rnor*step_delta
  END DO

  CALL LLdelta(deltanew, k, tau, lambda12, llnew)
  CALL LLdelta(delta, k, tau, lambda12, llold)

  lldiff = llnew - llold
  alpha = min(1.0_dp, exp(lldiff))

  IF (alpha .LT. 1.0_dp) THEN
    CALL random_number(unif)
  ELSE
    unif = 0.0_dp
  END IF

  IF (unif .LE. alpha) THEN
    delta = deltanew
    naccept = naccept + 1
  END IF
if (delta.LT. 0.20) delta=0.20

  ! write current draws to output arrays
  IF (mod(i1, keep) .EQ. 0) THEN
    betadraw((i1/keep),1:k) = beta
    lambda12draw((i1/keep),1:k) = lambda12
    deltadraw(i1/keep) = delta
    taudraw(i1/keep) = tau
  END IF

  rejrate = 1.0_dp - real(naccept,dp)/real(r,dp)

END DO

!=========================================================================

CONTAINS

!=========================================================================

! This code generates one draw from the standard normal 
! distribution. Note that more efficient code is possible
! when more than one normal draw is required.
! This code is based on the Box-Muller method.

! Output arguments:
!	- fn_val	: random draw from N(0,1) distribution

SUBROUTINE rnorm(fn_val)

IMPLICIT NONE

! Precision statement:
INTEGER, PARAMETER :: dp = KIND(1.0d0)

! Output arguments:
REAL(dp), INTENT(OUT) :: fn_val

! Internal arguments:
REAL(dp) :: pi
REAL(dp), DIMENSION(1:2) :: u

pi = 3.14159265358979323846_dp

CALL random_number(u)

fn_val = sqrt(-2*log(u(1))) * cos(2*pi*u(2))

END SUBROUTINE rnorm

!=========================================================================

! This code generates one random draw from the inverse Gaussian distribution.
! The algorithm is based on: Michael, Schucany & Haas (1976), Generating
! random variates using transformations with multiple roots, The
! American Statistician, 30(2), p. 88-90.

! Input arguments:
!	- mu		: mean parameter of the InvGaussian distribution
!	- lambda	: shape parameter of the InvGaussian distribution

! Output arguments:
!	- fn_val	: random InvGaussian variate

SUBROUTINE rInvGaussian (mu, lambda, fn_val)

IMPLICIT NONE

! Precision statement
INTEGER, PARAMETER :: dp = KIND(1.0d0)

! Input arguments:
REAL(dp), INTENT(IN) :: mu, lambda

! Output arguments:
REAL(dp), INTENT(OUT) :: fn_val

! Internal arguments:
REAL(dp) :: nu, q, z

CALL rnorm(nu)

nu = nu**2.0_dp
q = mu + (nu*mu**2.0_dp)/(lambda*2.0_dp) - &
    mu/(2.0_dp*lambda)*sqrt(4.0_dp*mu*lambda*nu + mu**2.0_dp*nu**2.0_dp)
nu = mu/(mu+q)

CALL random_number(z)

IF (z .LE. nu) THEN
    fn_val = q
ELSE
    fn_val = mu**2.0_dp/q
END IF

END SUBROUTINE rInvGaussian

!=========================================================================

! Implementation of the Lanczos approximation of the gamma
! function. Only use this code when the goal is to compute
! the LOGgamma, i.e. log(lancz_gamma(x)). Imprecise as 
! approximation for the gamma function for larger x.
! See:
! Lanczos, Cornelius (1964). "A Precision Approximation of the 
! Gamma Function". SIAM Journal on Numerical Analysis series B 
! (Society for Industrial and Applied Mathematics) 1: 86-96.

! Input arguments:
!   - x       : point to evaluate

! Output arguments:
!   - fn_val  : LOG of Lanczos approximation of the gamma function

RECURSIVE SUBROUTINE lancz_gamma(x, fn_val)

IMPLICIT NONE   

! Precision statement:
INTEGER, PARAMETER :: dp = KIND(1.0d0)

! Input arguments:
REAL(dp), INTENT(IN) :: x

! Output arguments:
REAL(dp), INTENT(OUT) :: fn_val

! Local arguments:
INTEGER :: i1
REAL(dp) :: t, w, a, b
REAL(dp), DIMENSION(1:8) :: c
INTEGER, PARAMETER :: cg = 7
REAL(dp), PARAMETER :: pi = 3.14159265358979324_dp
REAL(dp), DIMENSION(0:8), PARAMETER :: p = &
        (/ 0.99999999999980993_dp, 676.5203681218851_dp, -1259.1392167224028_dp, &
        771.32342877765313_dp, -176.61502916214059_dp, 12.507343278686905_dp, &
        -0.13857109526572012_dp, 9.9843695780195716D-6, 1.5056327351493116D-7 /)

a = x

IF (a < .5_dp) THEN
        CALL lancz_gamma(1.0_dp - a, b)
        fn_val = pi / (sin(pi*a) * b)
ELSE
        a = a - 1.0_dp
        c(1) = a + 1.0_dp
        DO i1 = 1,7
          c(i1+1) = c(i1) + 1.0_dp
        END DO
        t = p(0) + sum(p(1:8)/c)
        w = a + REAL(cg,dp) + .5_dp
        fn_val = sqrt(2.0_dp*pi) * w**(a+.5_dp) * exp(-w) * t
END IF

END SUBROUTINE lancz_gamma

!=========================================================================

! This code evaluates the logliklihood of the delta parameter

SUBROUTINE LLdelta (x, k, tau, lambda12, fn_val)

IMPLICIT NONE   

! Precision statement:
INTEGER, PARAMETER :: dp = KIND(1.0d0)

! Input arguments:
INTEGER, INTENT(IN) :: k
REAL(dp), INTENT(IN) :: x, tau
REAL(dp), INTENT(IN), DIMENSION(1:k) :: lambda12

! Output arguments:
REAL(dp), INTENT(OUT) :: fn_val

! Local arguments:
REAL(dp) :: w

CALL lancz_gamma(x, w)
fn_val = -real(k,dp)*log(w) + real(k,dp)*x*log(tau) - 2.0_dp*x* &
         sum(log(sqrt(lambda12)))

END SUBROUTINE LLdelta

!=========================================================================

! Generates one random draw from the gamma distribution with
! mean = shape*scale. The algorithm is based on Marsaglia & Tsang 
! "A Simple Method for Gererating Gamma Variables" (2000)

! Input arguments:
!	- shape		: shape parameter of the gamma distribution
!	- scale		: scale parameter of the gamma distribution

! Output arguments:
!	- fn_val	: random gamma variate Gamma(shape, scale)

SUBROUTINE rgamma (shape, scale, fn_val)

IMPLICIT NONE

! Precision statement:
INTEGER, PARAMETER :: dp = KIND(1.0d0)

! Input arguments:
REAL(dp), INTENT(IN) :: shape, scale

! Output arguments:
REAL(dp), INTENT(OUT) :: fn_val

! Internal arguments
REAL(dp) :: a, d, c, x, v, u
LOGICAL :: flag 


IF (shape < 1.0_dp) THEN
  a = shape + 1.0_dp
ELSE
  a = shape
END IF

d = a - 1.0_dp/3.0_dp
c = 1.0_dp/SQRT(9.0_dp*d)

flag = .TRUE.

DO WHILE (flag)
  v = 0.0_dp

  DO WHILE (v <= 0.0_dp)
    CALL rnorm(x)
    v = (1.0_dp + c*x)**3.0_dp
  END DO

  CALL RANDOM_NUMBER(u)

  IF (u < (1.0_dp-(0.0331_dp*(x**4.0_dp)))) THEN
    fn_val = d*v
    flag = .FALSE.
  END IF

  IF (LOG(u) < ((0.5_dp*x*x) + (d*(1.0_dp - v + LOG(v))))) THEN
    fn_val = d*v
    flag = .FALSE.
  END IF

END DO


IF (shape < 1.0_dp) THEN
  CALL RANDOM_NUMBER(u)
  fn_val = (fn_val * (u**(1.0_dp/shape))) * scale
ELSE
  fn_val = fn_val * scale
END IF


END SUBROUTINE rgamma

!=========================================================================


! This code evaluates the logliklihood of the lambdA2 parameter

SUBROUTINE lambd (x, k, nu12, nu2, beta, fn_val)

IMPLICIT NONE   

! Precision statement:
INTEGER, PARAMETER :: dp = KIND(1.0d0)

! Input arguments:
INTEGER, INTENT(IN) :: k
REAL(dp), INTENT(IN) :: x, nu12,nu2
REAL(dp), INTENT(IN), DIMENSION(1:k) :: beta

! Output arguments:
REAL(dp), INTENT(OUT) :: fn_val

! Local arguments:



fn_val =  -0.5_dp * x * nu12 * sum(beta**2.0_dp)/nu2

END SUBROUTINE lambd

!=========================================================================

END SUBROUTINE Aelasticnet
