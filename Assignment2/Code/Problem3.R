# A3, a
set.seed(1)
yA = c(12,9,12,14,13,13,15,8,15,6)
yB = c(11,11,10,9,9,8,7,10,6,8,8,9,7)
nA = 10
nB = 13
aA = 120
bA = 10
aB = 12
bB = 1
thetaA.posterior = rgamma(1000, aA+sum(yA), bA+nA)
thetaB.posterior = rgamma(1000, aB+sum(yB), bB+nB)
mean(thetaA.posterior > thetaB.posterior)

# A3, b
set.seed(2)
means = c()
n0s = 1:20
for (n0 in n0s)
{
aB = 12*n0
bB = 1*n0
thetaA.posterior1 = rgamma(10000, aA+sum(yA), bA+nA)
thetaB.posterior1 = rgamma(10000, aB+sum(yB), bB+nB)
avg = mean(thetaA.posterior1 > thetaB.posterior1)
means = append(means, avg)
}
plot(n0s, means, type="l", xlab="n0", ylab="Pr(θB < θA | yA,yB)")

# A3, c
set.seed(3)
yA = c(12,9,12,14,13,13,15,8,15,6)
yB = c(11,11,10,9,9,8,7,10,6,8,8,9,7)
nA = 10
nB = 13
aA = 120
bA = 10
aB = 12
bB = 1
thetaA.posterior = rgamma(10000, aA+sum(yA), bA+nA)
thetaB.posterior = rgamma(10000, aB+sum(yB), bB+nB)
yA.mc = rpois(10000,thetaA.posterior)
yB.mc = rpois(10000,thetaB.posterior)
mean(yA.mc > yB.mc)


set.seed(4)
yA = c(12,9,12,14,13,13,15,8,15,6)
yB = c(11,11,10,9,9,8,7,10,6,8,8,9,7)
nA = 10
nB = 13
aA = 120
bA = 10
means = c()
n0s = 1:20
for (n0 in n0s)
{
  aB = 12*n0
  bB = 1*n0
  thetaA.posterior2 = rgamma(10000, aA+sum(yA), bA+nA)
  thetaB.posterior2 = rgamma(10000, aB+sum(yB), bB+nB)
  yA.mc = rpois(10000,thetaA.posterior2)
  yB.mc = rpois(10000,thetaB.posterior2)
  avg = mean(yA.mc > yB.mc)
  means = append(means, avg)
}
plot(n0s, means, type="l", xlab="n0", ylab="Pr(YB < YA | yA,yB)")
