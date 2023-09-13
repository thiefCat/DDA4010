df = read.table("/Users/zhaosonglin/programs/R/DDA4010/Assignment/Assignment5/msparrownest.dat")
x = df$V2
y = df$V1
n = length(y) 
X = cbind(rep(1, n), x)
p = dim(X)[2]

rmvnorm <- function(n,mu,Sigma) {
  
  p<-length(mu)
  res<-matrix(0,nrow=n,ncol=p)
  if( n>0 & p>0 ) {
    E<-matrix(rnorm(n*p),n,p)
    res<-t(  t(E%*%chol(Sigma)) +c(mu))
  }
  return(res)
}


# (c)
# initialize
var.prop<- var(log(y+1/2))*solve(t(X)%*%X)
theta = rep(0,p)
Etheta.prior = c(0, 0)
SDtheta.prior = sqrt(c(25, 0.25))
S<-50000
THETA = matrix(NA,nrow=S,ncol=p)
ac<-0
set.seed(1)

## MCMC
for(s in 1:S) {
  theta.star = t(rmvnorm(1, theta, var.prop ))
  
  p = exp(X %*% theta.star) / (1+exp(X %*% theta.star))
  pp = exp(X %*% theta) / (1+exp(X %*% theta))
  r.log = sum(dbinom(y, 1, p, log=T)) + sum(dnorm(theta.star, Etheta.prior, SDtheta.prior, log = TRUE)) -
          sum(dbinom(y, 1, pp, log=T)) - sum(dnorm(theta, Etheta.prior, SDtheta.prior, log = TRUE))
  if( log(runif(1)) < r.log ) { 
    theta = theta.star 
    }
  THETA[s,] = theta
}

library(coda)
apply(THETA,2,effectiveSize)

# (d)
# alpha
xaxis1 <- seq(-15,15,length=1000)
plot(xaxis1, dnorm(xaxis1, 0, 5), type="l",col="red",xlab=expression(alpha), ylim=c(0,0.15), main="alpha", ylab="densidy")
lines(density(THETA[,1]), type="l",col="blue",xlab=expression(alpha), main="")
legend('topright',col=c("blue", "red"),legend=c('Posterior', 'Prior'), cex = 0.7,lty = 1)
# beta
xaxis2 <- seq(-2,2,length=100)
plot(xaxis2, dnorm(xaxis2, 0, 0.5), type="l",col="red",xlab=expression(alpha), ylim=c(0,2), main="beta", ylab="densidy")
lines(density(THETA[,2]), type="l",col="blue",xlab=expression(beta), main="")
legend('topright',col=c("blue", "red"),legend=c('Posterior', 'Prior'), cex = 0.7,lty = 1)


# (e)
l=50
x_s = seq(10, 15, length = l)
# M = matrix(NA, nrow=l,ncol=S)
QUANTILE = matrix(NA, nrow=l, ncol=3)
ALPHA = THETA[,1]
BETA = THETA[,2]
for (i in 1:l){
  xx = x_s[i]
  f = exp(ALPHA + BETA * xx) / ( 1+exp(ALPHA + BETA * xx) )
  QUANTILE[i,] = quantile(f, c(0.025, 0.5, 0.975))
}
plot(x_s, QUANTILE[, 2], type="l", ylim=c(0,1), xlab="wingspan",col="black")
lines(x_s, QUANTILE[, 1], type="l",col="blue")
lines(x_s, QUANTILE[, 3], type="l",col="red")
legend('topright',col=c("red", "black", "blue"),legend=c('97.5%', '50%', '2.5%'), cex = 0.6,lty = 1)
