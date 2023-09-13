df = read.table("/Users/zhaosonglin/programs/R/DDA4010/Assignment/Assignment4/interexp.dat", header=1)
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

# Question a
theta_A = mean(df$yA, na.rm = TRUE)
theta_B = mean(df$yB, na.rm = TRUE)
sigmasq_A = var(df$yA, na.rm = TRUE)
sigmasq_B = var(df$yB, na.rm = TRUE)
rou = cor(df$yA, df$yB, use = "complete.obs")

theta_A
theta_B
sigmasq_A
sigmasq_B
rou



# Question b
YA = df$yA
YB = df$yB
n = length(YA)

isNA = apply(df,2, function(x) is.na(x))
YA[isNA[,1]] = theta_A + (YB[isNA[,1]] - theta_B) * rou * sqrt(sigmasq_B/sigmasq_A)
YB[isNA[,2]] = theta_B + (YA[isNA[,2]] - theta_A) * rou * sqrt(sigmasq_A/sigmasq_B)
YA
YB

t.test(YA, YB, paired = TRUE, alternative = "two.sided")



# Question c
yA = df$yA
yB = df$yB
n = length(yA)
p = 2
S = 10000
THETA = matrix(0,S,p)
SIGMA = array(dim = c(p, p, S))
# initialize Y.full by YA, YB
# initialize Sigma by covariance of complete data
Sigma = cov(df, use = "complete.obs")
O = apply(df,2, function(x) is.na(x))
set.seed(1)
for(s in 1:S){
  Y = cbind(YA,YB)
  ybar = apply(Y, 2, mean)
  theta = rmvnorm(1,ybar,Sigma/(n+1))  
  
  S = ( t(Y)-c(ybar) ) %*% t( t(Y)-c(ybar) ) / n
  Sn = (n+1) * (( ( t(theta)-ybar ) %*% t( t(theta)-ybar ) ) + S ) 
  Sigma = solve( rwish(1, n+p+2, solve(Sn)) )
  
  for(i in 27:n)
  { 
    b = ( O[i,]==1 )
    a = ( O[i,]==0 )
    iSa =  solve(Sigma[a,a])
    beta.j = Sigma[b,a]%*%iSa
    s2.j = Sigma[b,b] - Sigma[b,a]%*%iSa%*%Sigma[a,b]
    theta.j = theta[b] + beta.j %*% (t(Y[i,a])-theta[a])
    Y[i,b] = rmvnorm(1,theta.j,s2.j)
  }
  YA = Y[,1]
  YB = Y[,2]
  THETA[s,] = theta
  SIGMA[,,s] = Sigma
}

plot(THETA[,1], THETA[,2],pch=".",xlab = "theta1", ylab = "theta2",main="scatter plot of posterior mean")

mean(THETA[,1]-THETA[,2])
quantile(THETA[,1]-THETA[,2], prob=c(.025,.975))

