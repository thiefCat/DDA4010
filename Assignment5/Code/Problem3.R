df = read.table("/Users/zhaosonglin/programs/R/DDA4010/Assignment/Assignment5/crime.dat", header=1)

rmvnorm <- function(n,mu,Sigma) {
  
  p<-length(mu)
  res<-matrix(0,nrow=n,ncol=p)
  if( n>0 & p>0 ) {
    E<-matrix(rnorm(n*p),n,p)
    res<-t(  t(E%*%chol(Sigma)) +c(mu))
  }
  return(res)
}

# bayesian estimates
X = data.matrix(df[,2:16])
y = df[,1]
n = dim(X)[1]
p = dim(X)[2]
S = 5000
g = n 
nu0 = 2
sigma0_sq = 1
SSRg = t(y) %*% (diag(1, n) - (g/(g + 1)) * X %*% solve(t(X) %*% X) %*% t(X)) %*% y
BETA = matrix(NA, nrow=S, ncol=p)
SIGMA_SQ = matrix(NA, nrow=S, ncol=1)
set.seed(1)
for (s in 1:S){
  sigma_sq = 1/rgamma(1, (nu0+n)/2, (nu0*sigma0_sq+SSRg)/2)
  Mbeta = (g/(g + 1)) * solve(t(X) %*% X) %*% t(X) %*% y
  Vbeta = (g/(g + 1)) * sigma_sq * solve(t(X) %*% X)
  beta = rmvnorm(1, Mbeta, Vbeta)
  BETA[s,] = beta
  SIGMA_SQ[s] = sigma_sq
}

# posterior mean
beta.posterior = apply(BETA, 2, mean)
beta.posterior
# posterior CI
apply(BETA, 2, quantile, probs = c(0.025, 0.975))

# least square estimates
beta.ols = solve(t(X)%*%X)%*%t(X)%*%y
beta.ols



# (b)
# (1)
set.seed(3)
indexs = sample.int(length(y), 24, replace=FALSE)
Xtr = X[indexs,]
Xte = X[-indexs,]
ytr = y[indexs]
yte = y[-indexs]
beta.ols = solve(t(Xtr)%*%Xtr)%*%t(Xtr)%*%ytr
y_predict = Xte %*% beta.ols
plot(yte, y_predict, ylab=expression(widehat(y)), xlab=expression(y), main="OLS")
abline(a = 0, b = 1)

error = sum((y_predict-yte)^2) / length(yte)
error



# (2)
n = dim(Xtr)[1]
p = dim(Xtr)[2]
S = 5000
g = n 
nu0 = 2
sigma0_sq = 1
SSRg = t(ytr) %*% (diag(1, n) - (g/(g + 1)) * Xtr %*% solve(t(Xtr) %*% Xtr) %*% t(Xtr)) %*% ytr
BETA = matrix(NA, nrow=S, ncol=p)
SIGMA_SQ = matrix(NA, nrow=S, ncol=1)
set.seed(1)
for (s in 1:S){
  sigma_sq = 1/rgamma(1, (nu0+n)/2, (nu0*sigma0_sq+SSRg)/2)
  Mbeta = (g/(g + 1)) * solve(t(Xtr) %*% Xtr) %*% t(Xtr) %*% ytr
  Vbeta = (g/(g + 1)) * sigma_sq * solve(t(Xtr) %*% Xtr)
  beta = rmvnorm(1, Mbeta, Vbeta)
  BETA[s,] = beta
  SIGMA_SQ[s] = sigma_sq
}

# posterior mean
beta.posterior = apply(BETA, 2, mean)
y_predict = Xte %*% beta.posterior
plot(yte, y_predict, ylab=expression(widehat(y)), xlab=expression(y), main="bayes")
abline(a = 0, b = 1)
error = sum((y_predict-yte)^2) / length(yte)
error



# (c)
set.seed(5)
Error_OLS = matrix(50,1)
Error_Bayes = matrix(50,1)
for (i in 1:50){
  # OLS
  indexs = sample.int(length(y), 24, replace=FALSE)
  Xtr = X[indexs,]
  Xte = X[-indexs,]
  ytr = y[indexs]
  yte = y[-indexs]
  beta.ols = solve(t(Xtr)%*%Xtr)%*%t(Xtr)%*%ytr
  y_predict = Xte %*% beta.ols
  error_ols = sum((y_predict-yte)^2) / length(yte)
  Error_OLS[i] = error_ols
  
  # Bayes
  n = dim(Xtr)[1]
  p = dim(Xtr)[2]
  S = 1000
  g = n 
  nu0 = 2
  sigma0_sq = 1
  SSRg = t(ytr) %*% (diag(1, n) - (g/(g + 1)) * Xtr %*% solve(t(Xtr) %*% Xtr) %*% t(Xtr)) %*% ytr
  BETA = matrix(NA, nrow=S, ncol=p)
  SIGMA_SQ = matrix(NA, nrow=S, ncol=1)
  for (s in 1:S){
    sigma_sq = 1/rgamma(1, (nu0+n)/2, (nu0*sigma0_sq+SSRg)/2)
    Mbeta = (g/(g + 1)) * solve(t(Xtr) %*% Xtr) %*% t(Xtr) %*% ytr
    Vbeta = (g/(g + 1)) * sigma_sq * solve(t(Xtr) %*% Xtr)
    beta = rmvnorm(1, Mbeta, Vbeta)
    BETA[s,] = beta
    SIGMA_SQ[s] = sigma_sq
  }
  
  # posterior mean
  beta.posterior = apply(BETA, 2, mean)
  y_predict = Xte %*% beta.posterior
  error_bayes = sum((y_predict-yte)^2) / length(yte)
  Error_Bayes[i] = error_bayes
}
mean(Error_OLS)
mean(Error_Bayes)
mean(Error_Bayes < Error_OLS)



