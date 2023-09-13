# (a)
set.seed(1)
Y <- scan("/Users/zhaosonglin/programs/R/DDA4010/Assignment/Assignment3/glucose.dat")
hist(Y, xlab = "concentration", main="distribution of the data")


# (c)
a = 1
b = 1
u0 = 120
tao0sq = 200
sigma0sq = 1000
v0 = 10

# initial values
n = length(Y)
p = 0.5
theta1 = mean(Y)
theta2 = mean(Y)
sigma1sq = var(Y)
sigma2sq = var(Y)
S = 10000
THETA_min = matrix(NA, S, 1)
THETA_max = matrix(NA, S, 1)
X_PRED = matrix(NA, S, 1)
Y_PRED = matrix(NA, S, 1)
for (s in 1:S){
  # xi
  p1 = p * dnorm(Y, theta1, sqrt(sigma1sq))
  p2 = (1-p) * dnorm(Y, theta2, sqrt(sigma2sq))
  prob = p1/(p1+p2)
  X = rbinom(n, 1, prob)  # X=1: prob, X=0: 1-Prob
  # p
  n1 = sum(X)
  n2 = n - n1
  p = rbeta(1, a+n1, b+n2)
  # theta1
  Y1 = Y[X == 1]
  Y2 = Y[X == 0]  # X==2
  y1bar = mean(Y1)
  y2bar = mean(Y2)
  tao1nsq = 1/(1/tao0sq + n1/sigma1sq)
  u1n = tao1nsq * (u0/tao0sq + n1*y1bar/sigma1sq)
  theta1 = rnorm(1, u1n, sqrt(tao1nsq))
  # theta2
  tao2nsq = 1/(1/tao0sq + n2/sigma2sq)
  u2n = tao2nsq * (u0/tao0sq + n2*y2bar/sigma2sq)
  theta2 = rnorm(1, u2n, sqrt(tao2nsq))
  # sigma1sq
  a1 = (v0+n1) / 2
  b1 = ( v0*sigma0sq+ sum( (Y1-theta1)^2 ) )/2
  sigma1sq = 1 / rgamma(1, a1, b1)
  # sigma2sq
  a2 = (v0+n2) / 2
  b2 = ( v0*sigma0sq+ sum( (Y2-theta2)^2 ) )/2
  sigma2sq = 1 / rgamma(1, a2, b2)
  THETA_min[s] = min(theta1, theta2)
  THETA_max[s] = max(theta1, theta2)
  
  x_pred = rbinom(1, 1, p)
  if (x_pred == 1){
    y_pred = rnorm(1, theta1, sqrt(sigma1sq))
  }
  if (x_pred == 0){
    y_pred = rnorm(1, theta2, sqrt(sigma2sq))
  }
  X_PRED[s] = x_pred
  Y_PRED[s] = y_pred

}

library(coda)

acf(THETA_min)
effectiveSize(THETA_min)
acf(THETA_max)
effectiveSize(THETA_max)


# (d)
# hist(Y_PRED, xlab = "concentration", main="empirical distribution of Y_predict and Y")
plot(density(Y_PRED), xlab="y", main="empirical distribution of Y_predict and Y", col="red")
lines(density(Y), col="blue")

