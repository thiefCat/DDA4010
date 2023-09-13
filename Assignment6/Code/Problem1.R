df = read.table("/Users/zhaosonglin/programs/R/DDA4010/Assignment/Assignment6/azdiabetes.dat", header=1)
set.seed(1)
y = (df[,8] == "Yes") * 1
X = data.matrix(df[,c(1,3,5,6,7)])
X = scale(X)
# X = cbind(rep(1, n), x)


rmvnorm <- function(n,mu,Sigma) {

  p<-length(mu)
  res<-matrix(0,nrow=n,ncol=p)
  if( n>0 & p>0 ) {
    E<-matrix(rnorm(n*p),n,p)
    res<-t(  t(E%*%chol(Sigma)) +c(mu))
  }
  return(res)
}

expit = function(x) {
  exp(x)/(1+exp(x))
}

fix = function(x) {
  ind = is.nan(x)|is.na(x) 
  x[ind] = 0
  x
}


p = 5
S = 10000
delta1 = 0.3
delta2 = 2

acs1 = 0
acs2 = 0

gamma = rbinom(p, 1, 0.5)
beta = rep(0, p)
beta0 = 0
BETA0 = matrix(NA,nrow=S,ncol=1)
BETA = matrix(NA,nrow=S,ncol=p)
GAMMA = matrix(NA,nrow=S,ncol=p)
for(s in 1:S) {
  # beta0
  beta0.star = rnorm(1, beta0, delta1)
  a1 = expit(as.matrix(X[,gamma==1]) %*% beta[gamma==1] + beta0.star)
  b1 = 1-a1
  a2 = expit(as.matrix(X[,gamma==1]) %*% beta[gamma==1] + beta0)
  b2 = 1-a2
  r.log = sum(y*log(a1) + (1-y)*log(b1)) + dnorm(beta0.star, 0, 4, log=T) -
    sum(y*log(a2) + (1-y)*log(b2)) - dnorm(beta0, 0, 4, log=T)
  
  # log(runif(1)) < r.log
  if( log(runif(1)) < r.log ) {
    beta0 = beta0.star
    acs1 = acs1+1
  }
  BETA0[s] = beta0

  # beta
  for(j in 1:p) {
    beta.star = beta
    beta.star[j] = rnorm(1, beta[j], delta2) 
    pi = expit(beta0 + as.matrix(X[,gamma==1]) %*% beta[gamma==1]) 
    pi_p = expit(beta0 + as.matrix(X[,gamma==1]) %*% beta.star[gamma==1]) 
    log_a = sum(fix(y*log(pi_p))+fix((1-y)*log(1-pi_p))) + dnorm(beta.star[j], 0, 2, log=T)
    log_b = sum(fix(y*log(pi))+fix((1-y)*log(1-pi))) + dnorm(beta[j], 0, 2, log=T)
    r.log = log_a - log_b
    if( log(runif(1)) < r.log ) {
      beta[j] = beta.star[j]
      acs2 = acs2+1
    }
  }
  BETA[s,] = beta

  # gamma
  for(j in sample(p)) {
    gamma_a = gamma_b = gamma 
    gamma_a[j] = 1 
    gamma_b[j] = 0 
    pi_a = expit(beta0 + as.matrix(X[,gamma_a==1]) %*% beta[gamma_a==1]) 
    log_a = sum(fix(y*log(pi_a))+fix((1-y)*log(1-pi_a)))
    pi_b = expit(beta0 + as.matrix(X[,gamma_b==1]) %*% beta[gamma_b==1]) 
    log_b = sum(fix(y*log(pi_b))+fix((1-y)*log(1-pi_b)))
    log_odds = log_a - log_b
    u = runif(1)
    gamma[j] = ifelse(u < expit(log_odds), 1, 0)
  }
  GAMMA[s,] = gamma
}
BETA0 = BETA0[5000:10000,]
BETA = BETA[5000:10000,]
GAMMA = GAMMA[5000:10000,]
ACR1 = acs1/S
ACR2 = acs2/(S*p)
ACR1
ACR2

par(mfrow = c(3,2))
plot(BETA0, type = 'l', main = 'beta0') 
plot(BETA[,1], type = 'l', main = 'beta1') 
plot(BETA[,2], type = 'l', main = 'beta2') 
plot(BETA[,3], type = 'l', main = ' beta3') 
plot(BETA[,4], type = 'l', main = 'beta4') 
plot(BETA[,5], type = 'l', main = 'beta5')

par(mfrow = c(3,2))
plot(BETA0, type = 'l', main = 'beta0') 
plot((BETA*GAMMA)[,1], type = 'l', main = 'beta1*gamma1') 
plot((BETA*GAMMA)[,2], type = 'l', main = 'beta2*gamma2') 
plot((BETA*GAMMA)[,3], type = 'l', main = ' beta3*gamma3') 
plot((BETA*GAMMA)[,4], type = 'l', main = 'beta4*gamma4') 
plot((BETA*GAMMA)[,5], type = 'l', main = 'beta5*gamma5')

library(coda)
apply(BETA, 2, effectiveSize)
apply(BETA*GAMMA, 2, effectiveSize)





# (b)
table(as.data.frame(GAMMA))
sort(table(as.data.frame(GAMMA)), decreasing = TRUE)[1:5] / 5000


# (c)
BETAGAMMA = (BETA*GAMMA)
plot(density(BETAGAMMA[,1]))
plot(density(BETAGAMMA[,2]))
plot(density(BETAGAMMA[,3]))
plot(density(BETAGAMMA[,4]))
plot(density(BETAGAMMA[,5]))

apply(BETAGAMMA, 2, mean)
apply(GAMMA, 2, mean)

