blue <- scan("/Users/zhaosonglin/programs/R/DDA4010/Assignment/Assignment4/bluecrab.dat")
orange <- scan("/Users/zhaosonglin/programs/R/DDA4010/Assignment/Assignment4/orangecrab.dat")
X = matrix(blue, 50, 2, byrow=1)
Y = matrix(orange, 50, 2, byrow=1)


rmvnorm <- function(n,mu,Sigma) {
  p = length(mu)
  if( n>0 & p>0 ) {
    E = matrix(rnorm(n*p),n,p)
    res = t(  t(E%*%chol(Sigma)) + c(mu))  
  }
  return(res)
}

rwish<-function(n,nu0,S0)
{
  sS0 = chol(S0)
  S = array( dim=c( dim(S0),n ) ) # p * p * n
  for(i in 1:n)
  {
    U = matrix(rnorm(nu0 * dim(S0)[1]), nu0, dim(S0)[1])
    Z = U %*% sS0
    S[,,i] = t(Z)%*%Z
  }
  return(S[,,1:n])
}

xbar = apply(X,2,mean)
ybar = apply(Y,2,mean)
x_cov = cov(X)
y_cov = cov(Y)

# prior
x_mu0 = xbar
x_L0 = x_cov
x_nu0 = 4
x_S0 = x_cov

y_mu0 = ybar
y_L0 = y_cov
y_nu0 = 4
y_S0 = y_cov
x_n = nrow(X)
y_n = nrow(Y)
p = ncol(X)

# gibbs sampling for bluecrab
set.seed(1)
S = 10000
x_THETA = matrix(0,S,p)
x_SIGMA = array(dim = c(p, p, S))
x_Sigma = x_cov
for(s in 1:S){
  
  x_Ln = solve( solve(x_L0) + x_n*solve(x_Sigma) )
  x_mun = x_Ln %*% ( solve(x_L0)%*%x_mu0 + x_n*solve(x_Sigma)%*%xbar )
  x_theta = rmvnorm(1,x_mun,x_Ln)  
  
  x_Sn =  x_S0 + ( t(X)-c(x_theta) ) %*% t( t(X)-c(x_theta) ) 
  Sigma = solve( rwish(1, x_nu0+x_n, solve(x_Sn)) )

  x_THETA[s,] = x_theta
  x_SIGMA[,,s] = Sigma
}

x_posterior_mean = apply(x_THETA, 2, mean)
plot(x_THETA[,1], x_THETA[,2],pch=".",xlab = "theta1", ylab = "theta2",main="scatter plot of posterior mean for blue")



# gibbs sampling for orangecrab
y_THETA = matrix(0,S,p)
y_SIGMA = array(dim = c(p, p, S))
y_Sigma = y_cov
for(s in 1:S){
  
  y_Ln = solve( solve(y_L0) + y_n*solve(y_Sigma) )
  y_mun = y_Ln %*% ( solve(y_L0)%*%y_mu0 + y_n*solve(y_Sigma)%*%ybar )
  y_theta = rmvnorm(1,y_mun,y_Ln)  
  
  y_Sn =  y_S0 + ( t(Y)-c(y_theta) ) %*% t( t(Y)-c(y_theta) ) 
  Sigma = solve( rwish(1, y_nu0+y_n, solve(y_Sn)) )
  
  y_THETA[s,] = y_theta
  y_SIGMA[,,s] = Sigma
}
y_posterior_mean = apply(y_THETA, 2, mean)
plot(y_THETA[,1], y_THETA[,2],pch=".",xlab = "theta1", ylab = "theta2", main="scatter plot of posterior mean for orange")
# points(y_posterior_mean, pch="19", col="orange", lwd="20")
x_posterior_mean
y_posterior_mean
quantile(x_THETA[,1]-y_THETA[,1], prob=c(.025,.975))
quantile(x_THETA[,2]-y_THETA[,2], prob=c(.025,.975))


# Question c
# blue crab
x_rou = c()
for(s in 1:S){
  covariance = x_SIGMA[,,s]
  rou = covariance[1,2] / (sqrt(covariance[1])*sqrt(covariance[4]))
  x_rou = c(x_rou, rou)
}

y_rou = c()
for(s in 1:S){
  covariance = y_SIGMA[,,s]
  rou = covariance[1,2] / (sqrt(covariance[1])*sqrt(covariance[4]))
  y_rou = c(y_rou, rou)
}

dx = density(x_rou)
plot(dx, col="blue", lwd=2,ylim=range(0:150), main="posterior density of correlation coefficient"
     , xlab="rou", ylab="density")

dy = density(y_rou)
lines(dy, col="orange",lwd=2)

mean(y_rou > x_rou)


