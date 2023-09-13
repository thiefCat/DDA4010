# A4, a
set.seed(1)
# Please change the directory!!!!!!
data_nobach <- scan("/Users/zhaosonglin/programs/R/DDA4010/Assignment/Assignment2/menchild30nobach.dat")
data_bach <- scan("/Users/zhaosonglin/programs/R/DDA4010/Assignment/Assignment2/menchild30bach.dat")
nA = length(data_bach)
nB = length(data_nobach)
theta_A = rgamma(5000, 2+sum(data_bach), 1+nA)
theta_B = rgamma(5000, 2+sum(data_nobach), 1+nB)
YA_predict = rpois(5000, theta_A)
YB_predict = rpois(5000, theta_B)

# posterior predictive distribution of YA
dist.A = (table(c(YA_predict,0:7))-1)/sum(table(YA_predict))
plot(dist.A, lwd=5, col="black", main="posterior predictive distribution of YA", ylab="pdf")
# posterior predictive distribution of YB
dist.B = (table(c(YB_predict,0:7))-1)/sum(table(YB_predict))
plot(dist.B, lwd=5, col="black", main="posterior predictive distribution of YB", ylab="pdf")


# A4, c
quantile(theta_B-theta_A, c(0.025, 0.975))
quantile(YB_predict-YA_predict, c(0.025, 0.975))

# A4, d
set.seed(2)
emdB = (table(c(data_nobach,0:7))-1)/sum(table(data_nobach))
plot(0:7-.1, emdB, lwd=5, col="gray", ylab="pdf", type="h", xlab="Y")
test = rpois(5000, 1.4)
emdtest = (table(c(test,0:7))-1)/sum(table(test))
points(0:7+.1, emdtest,lwd=5,col="black",type="h")
legend(3.5,.35,
       legend=c("Poisson(1.4)","empirical distribution"),
       lwd=c(2,2),col=c("black","gray"),bty="n",cex=1)

# A4, e
m = matrix(NA, 5000, 2)
for (s in 1:5000)
{
  theta = theta_B[s]
  samples = rpois(218, theta)
  m[s,1] = sum(samples == 0)
  m[s,2] = sum(samples == 1)
}
num_of_0 = m[,1]
num_of_1 = m[,2]
plot(num_of_0, num_of_1, lwd=1, pch = 19, cex=.3)
x = sum(data_nobach == 0)
y = sum(data_nobach == 1)
points(x,y, col="blue", pch = 19, cex = 1)

