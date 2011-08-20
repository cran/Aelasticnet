enetbqr <-
function(formula, tau=0.5, burnin = 5000, mcmc=25000,keep=1 ){

p=tau

r=mcmc
sig_shape=0.1
sig_rate=0.1
step_delta=0.1
x <- as.matrix(x)

# A model with an intercept
if(all(x[, 1] == 1)) (x=x) else(x=cbind(1,x))

# A model without an intercept  
if(formula == (y ~0+x) ) x=x[,-1]
if(formula == (y ~x+0) ) x=x[,-1]
if(formula == (y ~x-1) ) x=x[,-1]
if(formula == (y ~-1+x)) x=x[,-1]

x <- as.matrix(x)  
n <- nrow(x)
k <- ncol(x)

## checks
if (tau<=0 || tau>=1) stop ("invalid tau:  tau should be >= 0 and <= 1. \nPlease respecify tau and call again.\n")
if(!(all(is.finite(y)) || all(is.finite(x)))) stop("All values must be finite and non-missing")
if(n != length(y)) stop("length(y) and nrow(x) must be the same. \nPlease respecify the length of y and call again.\n")
if(!(all(is.finite(y)) || all(is.finite(x)))) stop("All values must be finite and non-missing")
if(all(y!=as.numeric(y>0))) stop ("y must be binary variable. \nPlease respecify y as y <- as.numeric(y>0) and call enetbqr(y~x) again. If not, please use the function enetqr(y~x).\n")



# Start of algorithm

    ## Assign correct variable types
    n <- as.integer(n)
    k <- as.integer(k)
    r <- as.integer(r)
    keep <- as.integer(keep)
     y <- as.double(y)
    tau <- as.double(tau)
    step_delta <- as.double(step_delta)
    x <- as.double(x)
    sig_shape <- as.double(sig_shape)
    sig_rate <- as.double(sig_rate)
    betadraw <- double(k*r/keep)
    lambda12draw <- double(k*r/keep)
    deltadraw <- double(r/keep)
    taudraw <- double(r/keep)
    sdraw <- double(k*r/keep)
    rejrate <- double(1)


    ## Call Fortran routine
    fn_val <- .Fortran("QRbAL", n, k, r, keep, y, p, step_delta, x, 
                        sig_shape, sig_rate, betadraw, lambda12draw, 
                        deltadraw, taudraw, rejrate,sdraw)

fitted=betadraw=matrix(fn_val[[11]], nrow=r/keep, ncol=k)[burnin:r,]
Top.model=NULL
if(k>=2){
for(j in 1:k){if (median(fitted[, j])==0 | quantile(fitted[, j],0.25)==0) fitted[,j]=0}
for(j in 1:k){if (median(fitted[, j])==0 | quantile(fitted[, j],0.75)==0) fitted[,j]=0}
for(j in 1:k){Top.model[j]=mean(fitted[,j])}
} else{Top.model=mean(fitted)}
 

    return(list(betadraw=matrix(fn_val[[11]], nrow=r/keep, ncol=k)[burnin:r,],
lambda12draw=matrix(fn_val[[12]], nrow=r/keep, ncol=k)[burnin:r,],sdraw=matrix(fn_val[[16]], nrow=r/keep, ncol=k)[burnin:r,]
            ))
}

