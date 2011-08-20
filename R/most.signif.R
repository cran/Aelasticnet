most.signif <-
function(sv,tau=0.5,burnin = 3000, mcmc = 10000){
if(all(y!=as.numeric(y>0)))stop("y must be binary variable. \nPlease  use the function model.selection(y~x,tau).\n")

MEANs=NULL                      #   posterior mean of sdraw (varianse beta in lasso penalty)
MEANb=NULL                      #   posterior mean of beta
MEANl=NULL                      #   posterior mean of lambda12

p=tau


fit <- enetbqr(y~0+x,tau=p, burnin, mcmc, keep=1)

fitted.s=as.matrix(fit$s)       #   variance matrix
fitted.b=as.matrix(fit$beta)    #   prameter estimate matrix
fitted.l=as.matrix(fit$l)       #   lambda12 matrix


  
n1=dim(fitted.b)[1]
k =dim(fitted.b)[2]


#  posterior variance
if(tau<=0.50){for(j in 1:k){MEANs[j]=1/mean(fitted.s[,j])} } 
             else {for(j in 1:k){MEANs[j]=mean(fitted.s[,j])} } 

#  posterior mean
for(j in 1:k){MEANb[j]=mean(fitted.b[,j])}   

cat("The most significant variables at tau=",tau,"\n" )

Most.sign=which(MEANs !=0)
Most.sign=order(Most.sign)

x.r=c(1:k)

mat=cbind(x.r,MEANs)
mat=cbind(mat,MEANb)


results<-cbind(mat[,1],mat[,2],mat[,3])
results=data.frame(R1=results[,1], R2=results[,2], R3=results[,3] )

Order.v=results[with(results, order(-R2 )), ]

Variable=Order.v[,1]
Poster.Mean=Order.v[,3]
df = data.frame(Variable,   Poster.Mean)

print(df[1:sv,])
}
