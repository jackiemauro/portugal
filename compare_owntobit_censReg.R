rm(list = ls())

# make data
a    <- 2    # structural parameter of interest
b    <- 1    # strength of instrument
rho  <- 0.5  # degree of endogeneity

N    <- 1000
z    <- rnorm(N)
res1 <- rnorm(N)
res2 <- res1*rho + sqrt(1-rho*rho)*rnorm(N)
x    <- z*b + res1
ys   <- x*a + res2
d    <- (ys>0) #dummy variable
y    <- d*ys

data = data.frame(x1 = rep(1,length(y)), x2 = x, y1 = y, z = z)
attach(data)

# own optimization
source("loglik_ivtobit.R")
start = list(sig_v = sd(res1), sig_u = sd(res2),
               tau1 = 0, beta1 = 0,beta2 = 0,
               pi1 = 0, pi2 = 0)
solved_tob = optim(unlist(start), loglik_tobitiv_norho, hessian=TRUE, control = list(maxit = 5000))
solved_tob$par

# premade optimization
require(censReg)
reduced.form <- lm(x~z)
consistent.tobit <- censReg(y~fitted(reduced.form)+residuals(reduced.form))
summary(consistent.tobit)

#### with exogenous variable (not intercept) ####
# works fine if you don't include an x1 of all 1's and keep using "norho" fn

rm(list = ls())
detach(data)
a    <- 2    # structural parameter of interest
b    <- 1    # strength of instrument
rho  <- 0.5  # degree of endogeneity
N    <- 1000
z    <- rnorm(N)
res1 <- rnorm(N)
res2 <- res1*rho + sqrt(1-rho*rho)*rnorm(N)
b1   <- -1
a1   <- .5
exog <- 2 + .5*rnorm(N)
x2  <- z*b + exog*b1 + res1
ys   <- x2*a + exog*a1 + res2
d    <- (ys>0) #dummy variable
y1   <- d*ys

x1 = t(as.matrix(exog, ncol = 1))

# own optimization
source("loglik_ivtobit.R")
start = list(sig_v = sd(res1), sig_u = sd(res2),
             tau1 = 0, beta1 = 1,beta2 = 0,
             pi1 = 1, pi2 = 0)
solved_tob = optim(unlist(start), loglik_tobitiv_norho, hessian=TRUE, control = list(maxit = 5000))
solved_tob$par

# premade optimization
require(censReg)
reduced.form <- lm(x2~z+exog)
consistent.tobit <- censReg(y1~fitted(reduced.form)+exog+residuals(reduced.form))
summary(consistent.tobit)

#### with mult exogenous variables ####
# does fine if you don't include intercept. In this case, dgp has no intercept
rm(list = ls())
detach(data)
a    <- 2    # structural parameter of interest
b    <- 1    # strength of instrument
dg   <- 0.5  # degree of endogeneity
N    <- 1000
z    <- rnorm(N)
res1 <- rnorm(N)
res2 <- res1*dg + sqrt(1-dg*dg)*rnorm(N)
b11  <- -1
b12  <- 2
a11  <- .5
a12  <- .6
ex1  <- 2 + .5*rnorm(N)
ex2  <- 1 + rnorm(N)
x2   <- z*b + ex1*b11 + ex2*b12 + res1
ys   <- x2*a + ex1*a11 + ex2*a12 + res2
d    <- (ys>0) #dummy variable
y1   <- d*ys

# with intercept, does a bad job
x1 = rbind(rep(1,length(y1)), ex1, ex2)

# without intercept, does well
x1 = rbind(ex1,ex2)

# own optimization
source("loglik_ivtobit.R")
start = list(sig_v = sd(res1), sig_u = sd(res2),
             beta1 = c(2,2),beta2 = 2,
             pi1 = c(1,1), pi2 = 1)
solved_tob = optim(unlist(start), loglik_tobitiv_mat_norho, hessian=TRUE, control = list(maxit = 5000))
solved_tob$par

# premade optimization
require(censReg)
reduced.form <- lm(x2~z+ex1+ex2)
consistent.tobit <- censReg(y1~fitted(reduced.form)+ex1+ex2+residuals(reduced.form))
summary(consistent.tobit)


#### with mult exogenous variables and intercept ####
# doing ok, not as well as built in function
rm(list = ls())

a    <- 2    # structural parameter of interest
b    <- 1    # strength of instrument
dg   <- 0.5  # degree of endogeneity
N    <- 5000
z    <- rnorm(N)
res1 <- rnorm(N)
res2 <- res1*dg + sqrt(1-dg*dg)*rnorm(N)

# pi terms
int1 <- 3
b11  <- -1
b12  <- 2

# beta terms
int2 <- -1
a11  <- .5
a12  <- .6

ex1  <- 2 + .5*rnorm(N)
ex2  <- 1 + rnorm(N)
x2   <- int1 + z*b + ex1*b11 + ex2*b12 + res1
ys   <- int2 + x2*a + ex1*a11 + ex2*a12 + res2
d    <- (ys>0) #dummy variable
y1   <- d*ys

x1 = rbind(rep(1,length(y1)), ex1, ex2)

# own optimization
#source("loglik_ivtobit.R")
start = list(sig_v = sd(res1), sig_u = sd(res2),
             beta1 = c(2,2,2),beta2 = 2,
             pi1 = c(1,1,1), pi2 = 1)
true = list(sig_v = 1, sig_u = 1,
             beta1 = c(-1,.5,.6),beta2 = 2,
             pi1 = c(3,-1,2), pi2 = 1)
solved_tob = optim(unlist(start), loglik_tobitiv_mat_norho, hessian=TRUE, control = list(maxit = 5000))
solved_tob$par

# premade optimization
require(censReg)
reduced.form <- lm(x2~z+ex1+ex2)
consistent.tobit <- censReg(y1~fitted(reduced.form)+ex1+ex2+residuals(reduced.form))
summary(consistent.tobit)