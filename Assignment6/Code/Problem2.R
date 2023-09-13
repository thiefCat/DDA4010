df = read.table("/Users/zhaosonglin/programs/R/DDA4010/Assignment/Assignment6/pdensity.dat", header=1)
y = df$yield
x = df$density


rmvnorm<-function(n,mu,Sigma)
{ 
  E = matrix(rnorm(n*length(mu)),n,length(mu))
  t(  t(E%*%chol(Sigma)) +c(mu))
}


rwish<-function(n,nu0,S0)
{
  sS0 = chol(S0)
  S<-array( dim=c( dim(S0),n ) )
  for(i in 1:n)
  {
    Z <- matrix(rnorm(nu0 * dim(S0)[1]), nu0, dim(S0)[1]) %*% sS0
    S[,,i]<- t(Z)%*%Z
  }
  S[,,1:n]
}


n = 8
m = 10
# j = 2

# OLS
BETA.OLS = matrix(NA, 10, 3)
SIGMA_SQ.OLS = matrix(NA, 10, 1)
Y = list() 
X = list()
for (j in (1:10)){
  # yj = y[(8*(j-1)+1):(8*j)]
  # xj = x[(8*(j-1)+1):(8*j)]
  # Xj = cbind(rep(1,n),xj,xj**2)
  # betaj.ols = solve(t(Xj)%*%Xj)%*%t(Xj)%*%yj
  # BETA.OLS[j,] = betaj.ols  
  # SIGMA_SQ.OLS[j] = var(yj)
  
  Y[[j]] = y[(8*(j-1)+1):(8*j)]
  xj = x[(8*(j-1)+1):(8*j)]
  X[[j]] = cbind(rep(1,n),xj,xj**2)
  betaj.ols = solve(t(X[[j]])%*%X[[j]])%*%t(X[[j]])%*%Y[[j]]
  BETA.OLS[j,] = betaj.ols  
  SIGMA_SQ.OLS[j] = var(Y[[j]])
  
}

# plot
x = seq(0,10,length=100)
y = BETA.OLS[j,1] + BETA.OLS[j,2]*x + BETA.OLS[j,3]*x^2
plot(x,y,type="l",ylim=range(1:10))
for (j in (2:10)){
  x = seq(0,10,length=100)
  y = BETA.OLS[j,1] + BETA.OLS[j,2]*x + BETA.OLS[j,3]*x^2
  lines(x,y,type="l")
}

apply(BETA.OLS, 2, mean)
cov(BETA.OLS)
apply(SIGMA_SQ.OLS, 2, mean)


# (b) (c)
p<-dim(X[[1]])[2]
mu0 = apply(BETA.OLS, 2, mean)
L0 = cov(BETA.OLS) 
nu0 = 2 
s20 = apply(SIGMA_SQ.OLS, 2, mean)
eta0 = 4
S0 = cov(BETA.OLS) 
# Sigma.ps<-matrix(0,p,p)
# SIGMA.PS<-NULL
# BETA.ps<-BETA*0
# BETA.pp<-NULL
set.seed(1)

# initialize
theta = mu0
Sigma = cov(BETA.OLS) 
s2 = apply(SIGMA_SQ.OLS, 2, mean)
BETA = BETA.OLS
S = 10000
THETA.p = matrix(NA, S, p)
BETA.p = array(NA, dim = c(10,3,S))
SIGMA.p = array(NA, dim = c(3,3,S))
S2.p = matrix(NA, S, 1)
iL0 = solve(L0) 
iSigma = solve(Sigma)

for(s in 1:S) {
  # beta
  for(j in 1:m) 
  {  
    bj = solve( iSigma + t(X[[j]])%*%X[[j]]/s2 )
    aj = bj%*%( iSigma%*%theta + t(X[[j]]) %*% Y[[j]]/s2 )
    BETA[j,]<-rmvnorm(1,aj,bj) 
  } 
  BETA.p[,,s] = BETA
  
  # theta
  b =  solve( iL0 +  m*iSigma )
  a = b%*%( iL0%*%mu0 + iSigma%*%(m*mu0))
  theta = t(rmvnorm(1,a,b))
  THETA.p[s,] = theta
  
  # Sigma
  mtheta<-matrix(theta,m,p,byrow=TRUE)
  Sigma = solve( rwish(1, eta0+m, solve( S0+t(BETA-mtheta)%*%(BETA-mtheta) ) ) )
  SIGMA.p[,,s] = Sigma
  
  # s2
  RSS<-0
  for(j in 1:m) { 
    RSS = RSS+sum( (Y[[j]]-X[[j]]%*%BETA[j,] )^2 ) 
    }
  s2 = 1/rgamma(1,(nu0+m*n)/2, (nu0*s20+RSS)/2 )
  S2.p[s] = s2
}

# MCMC diagnosis
plot(S2.p, type = 'l', main = 'beta0') 
effectiveSize(S2.p)


# posterior mean of beta
BETA.p.mean = apply(BETA.p, c(1,2), mean)
x = seq(0,10,length=100)
y = BETA.p.mean[1,1] + BETA.p.mean[1,2]*x + BETA.p.mean[1,3]*x^2
plot(x,y,type="l",ylim=range(1:10))
for (j in (2:10)){
  x = seq(0,10,length=100)
  y = BETA.p.mean[j,1] + BETA.p.mean[j,2]*x + BETA.p.mean[j,3]*x^2
  lines(x,y,type="l")
}



# (d)
par(mfrow = c(1,3))
plot(density(THETA.p[,1]), main="posterior of theta1")
plot(density(THETA.p[,2]), main="posterior of theta2")
plot(density(THETA.p[,3]), main="posterior of theta3")

par(mfrow = c(3,3))
for (i in (1:3)){
  for (j in (1:3)){
    plot(density(SIGMA.p[i,j,]), main="posterior of sigma")
  }
}


# (e)
avgtheta = apply(THETA.p, 2, mean)
# y_avg = avgtheta[1] + avgtheta[2]*x + avgtheta[3]*x^2
x_max = - avgtheta[2] / (2*avgtheta[3])
x_max

X_max = - THETA.p[,2] / (2*THETA.p[,3])
quantile(X_max, c(0.025, 0.975))


