school1 <- scan("/Users/zhaosonglin/programs/R/DDA4010/Assignment/Assignment5/school1.dat")
school2 <- scan("/Users/zhaosonglin/programs/R/DDA4010/Assignment/Assignment5/school2.dat")
school3 <- scan("/Users/zhaosonglin/programs/R/DDA4010/Assignment/Assignment5/school3.dat")
school4 <- scan("/Users/zhaosonglin/programs/R/DDA4010/Assignment/Assignment5/school4.dat")
school5 <- scan("/Users/zhaosonglin/programs/R/DDA4010/Assignment/Assignment5/school5.dat")
school6 <- scan("/Users/zhaosonglin/programs/R/DDA4010/Assignment/Assignment5/school6.dat")
school7 <- scan("/Users/zhaosonglin/programs/R/DDA4010/Assignment/Assignment5/school7.dat")
school8 <- scan("/Users/zhaosonglin/programs/R/DDA4010/Assignment/Assignment5/school8.dat")

# (a)
Y<-list()
Y[[1]] = school1
Y[[2]] = school2
Y[[3]] = school3
Y[[4]] = school4
Y[[5]] = school5
Y[[6]] = school6
Y[[7]] = school7
Y[[8]] = school8
n = ybar = sample_var = rep(0, 8)
for (i in 1:8){
  n[i] = length(Y[[i]])
  ybar[i] = mean(Y[[i]])
  sample_var[i] = var(Y[[i]])
}
# priors
mu0 = 7; gamma0_sq = 5
tau0_sq = 10; eta0 = 2
nu0 = 2; sigma0_sq = 15

# initial values
theta = ybar
mu = mean(ybar)
sigma_sq = mean(sample_var)
tau_sq = var(theta)
m = 8

# store posterior values
S = 5000
THETA = matrix( nrow=S,ncol=m)
TAU_SQ = matrix( nrow=S,ncol=1)
MU = matrix( nrow=S,ncol=1)
SIGMA_SQ = matrix( nrow=S,ncol=1)
set.seed(1)
for(s in 1:S) 
  {
  # theta
  for (j in 1:m)
    {
    b = 1/(n[j]/sigma_sq+1/tau_sq)
    a = b*(ybar[j]*n[j]/sigma_sq+mu/tau_sq)
    theta[j] = rnorm(1,a,sqrt(b))
  }
  # mu
  Vmu = 1/(m/tau_sq+1/gamma0_sq)
  Emu = Vmu*(m*mean(theta)/tau_sq + mu0/gamma0_sq)
  mu = rnorm(1,Emu,sqrt(Vmu)) 
  # tau_sq
  a = (eta0+m) / 2
  b = ( eta0*tau0_sq + sum( (theta-mu)^2 ) ) / 2
  tau_sq = 1/rgamma(1,a,b)
  # sigma_sq
  a = (nu0 + sum(n)) / 2
  ss = 0
  for (j in 1:m){
    ss = ss + sum( (Y[[j]]-theta[j])^2 )
  }
  b = (nu0*sigma0_sq + ss) / 2
  sigma_sq<-1/rgamma(1,a,b)
  
  # store the posterior samples
  THETA[s,] = theta
  TAU_SQ[s] = tau_sq
  MU[s] = mu
  SIGMA_SQ[s] = sigma_sq
}

# trace plot
par(mfrow = c(3,1))
plot(MU, type = "l", xlab="iteration",ylab=expression(mu))
plot(TAU_SQ, type = "l", xlab="iteration",ylab=expression(tau^2))
plot(SIGMA_SQ, type = "l", xlab="iteration",ylab=expression(sigma^2))

# stationary plot
stationarity.plot<-function(x,...){
  
  S<-length(x)
  scan<-1:S
  ng<-min( round(S/100),10)
  group<-S*ceiling( ng*scan/S) /ng
  
  boxplot(x~group,...)               }
par(mfrow = c(1, 3))
stationarity.plot(MU,xlab="iteration",ylab=expression(mu))
stationarity.plot(TAU_SQ,xlab="iteration",ylab=expression(tau^2))
stationarity.plot(SIGMA_SQ,xlab="iteration",ylab=expression(sigma^2))

# effective sample size
library(coda)
effectiveSize(MU)
effectiveSize(TAU_SQ)
effectiveSize(SIGMA_SQ)


# (b)
# posterior means:
M = cbind(MU,TAU_SQ,SIGMA_SQ)
apply(M, MARGIN=2, FUN=mean)

# posterior confidence regions
apply(M, MARGIN = 2, FUN = quantile, probs = c(0.025, 0.975))

# prior densities and posterior densities
dinvgamma<-function(x,a,b) {
  ld<- a*log(b) -lgamma(a) -(a+1)*log(x)  -b/x
  exp(ld)
}
par(mfrow = c(1, 3))
plot(density(MU), type="l",col="blue", xlab=expression(mu),main="")
x1 = seq(3,13,length=100)
lines(x1, dnorm(x1, mu0, sqrt(gamma0_sq)), type="l",col="red")

plot(density(TAU_SQ), type="l",col="blue", xlab=expression(tau^2),main="")
x2 = seq(0,60,length=100)
lines(x2, dinvgamma(x2, eta0/2, eta0*tau0_sq/2), type="l",col="red")

plot(density(SIGMA_SQ), type="l",col="blue",xlab=expression(sigma^2),main="")
x3 = seq(5,30,length=100)
lines(x3, dinvgamma(x3, nu0/2, nu0*sigma0_sq/2), type="l",col="red")
legend('topright',col=c("blue", "red"),legend=c('Posterior', 'Prior'), cex = 0.7,lty = 1)



# (c)
par(mfrow=c(1,1))
R.posterior = TAU_SQ/(SIGMA_SQ+TAU_SQ)
plot(density(R.posterior), type="l",col="blue",xlab=expression(tau^2/(sigma^2 + tau^2)),main="prior and posterior of R")
TAU0_SQ = 1/rgamma(S, eta0/2, eta0*tau0_sq/2) 
SIGMA0_SQ = 1/rgamma(S, nu0/2, nu0*sigma0_sq/2)
R.prior = TAU0_SQ/(SIGMA0_SQ+TAU0_SQ)
lines(density(R.prior), type="l",col="red")
legend('topright',col=c("blue", "red"),legend=c('Posterior', 'Prior'), cex = 0.7,lty = 1)

mean(R.posterior)

# (d)
mean(THETA[,7]<THETA[,6])
minimum = apply(THETA, MARGIN=1, FUN=which.min)
mean(minimum==7)


# (e)
plot(ybar, apply(THETA, MARGIN=2, FUN=mean), cex=1, xlab=expression(bar(y)),ylab=expression(widehat(theta)))
abline(a = 0, b = 1)

# sample mean
sum = 0
for (j in 1:8){
  sum = sum + sum(Y[[j]])
}
# global sample mean
sum/sum(n)
# posterior mean of u
mean(MU)

