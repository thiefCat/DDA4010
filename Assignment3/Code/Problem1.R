# (a)
set.seed(1)
school1 <- scan("/Users/zhaosonglin/programs/R/DDA4010/Assignment/Assignment3/school1.dat")
school2 <- scan("/Users/zhaosonglin/programs/R/DDA4010/Assignment/Assignment3/school2.dat")
school3 <- scan("/Users/zhaosonglin/programs/R/DDA4010/Assignment/Assignment3/school3.dat")

# school 1
u0 <- 5
k0 <- 1
sigmasq0 <- 4
v0 <- 2

n_1 <- length(school1)
ybar_1 <- mean(school1)
s2_1 <- var(school1)   # sample variance

kn_1 <- k0 + n_1
vn_1 <- v0 + n_1
un_1 <- (k0 * u0 + n_1 * ybar_1) / kn_1
sigmasqn_1 <- (v0 * sigmasq0 + (n_1-1) * s2_1 + k0 * n_1 / kn_1 * (ybar_1 - u0)^2) / vn_1
# posterior means
un_1
# 95% CI for mean and standard deviation
print("CI for mean:")
ybar_1 + sqrt(s2_1) / sqrt(n_1) * qt(c(0.025, 0.975), n_1-1)
print("CI for sd:")
sqrt(1 / qgamma(c(0.975, 0.025), vn_1/2, sigmasqn_1 * vn_1/2))

# school 2
n_2 <- length(school2)
ybar_2 <- mean(school2)
s2_2 <- var(school2)   # sample variance

kn_2 <- k0 + n_2
vn_2 <- v0 + n_2
un_2 <- (k0 * u0 + n_2 * ybar_2) / kn_2
sigmasqn_2 <- (v0 * sigmasq0 + (n_2-1) * s2_2 + k0 * n_2 / kn_2 * (ybar_2 - u0)^2) / vn_2
# posterior means
un_2
# 95% CI for mean and standard deviation
print("CI for mean:")
ybar_2 + sqrt(s2_2) / sqrt(n_2) * qt(c(0.025, 0.975), n_2-1)
print("CI for sd:")
sqrt(1 / qgamma(c(0.975, 0.025), vn_2/2, sigmasqn_2 * vn_2/2))

# school 3
n_3 <- length(school3)
ybar_3 <- mean(school3)
s2_3 <- var(school3)   # sample variance

kn_3 <- k0 + n_3
vn_3 <- v0 + n_3
un_3 <- (k0 * u0 + n_3 * ybar_3) / kn_3
sigmasqn_3 <- (v0 * sigmasq0 + (n_3-1) * s2_3 + k0 * n_3 / kn_3 * (ybar_3 - u0)^2) / vn_3
# posterior means
un_3
# 95% CI for mean and standard deviation
print("CI for mean:")
ybar_3 + sqrt(s2_3) / sqrt(n_3) * qt(c(0.025, 0.975), n_3-1)
print("CI for sd:")
sqrt(1 / qgamma(c(0.975, 0.025), vn_3/2, sigmasqn_3 * vn_3/2))




# (b)
set.seed(2)
S <- 100000
sigmasq_1.post <- 1/rgamma(S, vn_1/2, sigmasqn_1 * vn_1/2)
sigmasq_2.post <- 1/rgamma(S, vn_2/2, sigmasqn_2 * vn_2/2)
sigmasq_3.post <- 1/rgamma(S, vn_3/2, sigmasqn_3 * vn_3/2)
u_1.post <- rnorm(S, un_1, sqrt(sigmasq_1.post/kn_1))
u_2.post <- rnorm(S, un_2, sqrt(sigmasq_2.post/kn_2))
u_3.post <- rnorm(S, un_3, sqrt(sigmasq_3.post/kn_3))
mean(u_1.post < u_2.post & u_2.post < u_3.post)
mean(u_1.post < u_3.post & u_3.post < u_2.post)
mean(u_2.post < u_1.post & u_1.post < u_3.post)
mean(u_2.post < u_3.post & u_3.post < u_1.post)
mean(u_3.post < u_2.post & u_2.post < u_1.post)
mean(u_3.post < u_1.post & u_1.post < u_2.post)



# (c)
S <- 10000
Y1_pred = rnorm(S, u_1.post, sqrt(sigmasq_1.post))
Y2_pred = rnorm(S, u_2.post, sqrt(sigmasq_2.post))
Y3_pred = rnorm(S, u_3.post, sqrt(sigmasq_3.post))
mean(Y1_pred < Y2_pred & Y2_pred < Y3_pred)
mean(Y1_pred < Y3_pred & Y3_pred < Y2_pred)
mean(Y2_pred < Y1_pred & Y1_pred < Y3_pred)
mean(Y2_pred < Y3_pred & Y3_pred < Y1_pred)
mean(Y3_pred < Y2_pred & Y2_pred < Y1_pred)
mean(Y3_pred < Y1_pred & Y1_pred < Y2_pred)



# (d)
mean(u_1.post > u_2.post & u_1.post > u_3.post)
mean(Y1_pred > Y2_pred & Y1_pred > Y3_pred)


