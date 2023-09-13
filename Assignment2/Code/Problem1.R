a1=1
theta1=1
eq1 = function(y){2/gamma(a1)*theta1^(2*a1)*y^(2*a1-1)*exp(-theta1^2*y^2)}
a2=2
theta2=1
eq2 = function(y){2/gamma(a2)*theta2^(2*a2)*y^(2*a2-1)*exp(-theta2^2*y^2)}
a3=1
theta3=2
eq3 = function(y){2/gamma(a3)*theta3^(2*a3)*y^(2*a3-1)*exp(-theta3^2*y^2)}
x = seq(0, 10, by=0.1)
plot(x,eq1(x),type="l",xlim=c(0.1,10),ylim=c(0,2), col="black", lwd=2, ylab="y")
lines(x,eq2(x), type="l", col="blue", lwd=2)
lines(x,eq3(x), type="l", col="green", lwd=2)
legend(6,1.5, legend=c("Galenshore(1,1)","Galenshore(2,1)","Galenshore(1,2)")
,lwd=c(2,2,2),col=c("black","blue","green"),bty="n",cex=1)