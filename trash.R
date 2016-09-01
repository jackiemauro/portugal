#####
# trash/testing file
####

setwd("C:/Users/jackie/Dropbox/portugal/theory paper")

rm(list = ls())
par(cex = 1.25)
par.old <- par()

img = "C:/Users/jackie/Dropbox/portugal/theory paper/JM methods paper/"
tab = "C:/Users/jackie/Dropbox/portugal/theory paper/JM methods paper/"

source("JM methods paper/data_simulations_mult.R")
source("start_values.R")
source("loglik_lgnorm.R")
source("loglik_cragg.R")
source("loglik_ivtobit.R")
source("comparison_function.R")
source("parameter_values.R")

setwd("C:/Users/jackie/Dropbox/portugal/theory paper")
source("hurdleIVClass.R")
source("startValuesFn.R")

require(MASS)
require(truncnorm)
require(xtable)

# reg = hurdleIV(y1 ~ x11 + x21 + x22 + x23,
#              exog = list(x11,x12),
#              start_val = F,
#              endog = list(x21~x11+x12+z1+z2+z3,
#                           x22~x11+z1+z2+z3,
#                           x23~x12+z1+z2+z3))


gen_cragg1(n = 1000, rho = F, params = cr_params_mult)
y1 = c(y1)
x11 = c(x1[1,]); x12 = c(x1[2,])
x21 = c(x2[1,]); x22 = c(x2[2,]); x23 = c(x2[3,])
z1 = c(z[1,]); z2 = c(z[2,]); z3 = c(z[3,])

regs = list(formula = y1 ~ x11 + x21 + x22 + x23,
            exog = list(x11,x12),
            inst = list(z1,z2,z3),
            endog = list(x21,x22,x23),
            start_val = F,
            endogReg = list(x21~x11+x12+z1+z2+z3,
                         x22~x11+z1+z2+z3,
                         x23~x12+z1+z2+z3))

### get the start values once
t = hurdleIV.start_vals(regs)

### optimization will run through the rest a few times
# get covariance
params = t
Sig = hurdleIV.covMat(params)

# get the means
j = length(regs$endog)

# first get the x2 means so they're in the right shape
mu_x2 = matrix(c(rep(NA,j*length(y1))),ncol = 3)
for(i in 1:j){
  formula = regs$endogReg[[i]]
  mf = model.frame(formula = formula)
  m = dim(mf)[2]
  x <- model.matrix(attr(mf, "terms"), data=mf)
  x1temp = x[,1:(m-j)]; ztemp = x[,(m-j+1):m]
  mu_x2[,i] = x1temp%*%params$pi1[[i]] + ztemp%*%params$pi2[[i]]
}

# then get the means for the y regressions -- this gets messy: clean up and check
formula = regs$formula
mf = model.frame(formula = formula)
m = dim(mf)[2]
x <- model.matrix(attr(mf,"terms"),data = mf)
x1 = x[,1:(m-j)]; x2 = x[,(m-j+1):m]
mu_y0 = x1%*%params$gamma1 + mu_x2%*%params$gamma2
mu_y1 = x1%*%params$beta1 + mu_x2%*%params$beta2

# now get the conditional means
#Parameters for x2
sig2_x2 = Sig[(m-j+1):m,(m-j+1):m]

#Parameters for y0star given x2
mu_y0_x2 = mu_y0 + t(Sig[1,(m-j+1):m]%*%solve(sig2_x2)%*%t(x2-mu_x2))
sig2_y0_x2 = Sig[1,1] - Sig[1,(m-j+1):m]%*%solve(sig2_x2)

#Parameters for log(y1star) given x2
mu_y1_x2 = mu_y1 + t(Sig[2,(m-j+1):m]%*%solve(sig2_x2)%*%t(x2-mu_x2))
sig2_y1_x2 = Sig[2,2] - Sig[2,(m-j+1):m]%*%solve(sig2_x2)

#Parameters for y0star given y1star and x2
mu_y0_y1x2 = mu_y0 + t(Sig[1,2:m,drop=FALSE]%*%solve(Sig[2:m,2:m])%*%t(cbind(y1-mu_y1,x2-mu_x2)))
sig2_y0_y1x2 = Sig[1,1] - Sig[1,2:m,drop=FALSE]%*%solve(Sig[2:m,2:m])%*%Sig[2:m,1,drop=FALSE]

###Calculate the contributions to the log likelihood.

#When y1=0:
ll0 = pnorm(0,mean=mu_y0_x2,sd=sqrt(sig2_y0_x2),log.p = TRUE) + 
  dnorm(x2,mean=mu_x2,sd=sqrt(sig2_x2),log = TRUE)

#When y1>0:
ll1 = pnorm(0,mean=mu_y0_y1x2, sd=sqrt(sig2_y0_y1x2),log.p=TRUE, lower.tail = FALSE) +
  dnorm(logy1, mean=mu_y1_x2, sd=sqrt(sig2_y1_x2),log = TRUE) + 
  dnorm(x2,mean=mu_x2,sd=sqrt(sig2_x2),log = TRUE)

#Combine them, based on y1
ll = ifelse(censored,ll0,ll1)

#Return the negative log-likelihood, since I will be using optim (which minimizes instead of maximizing)
-sum(ll)
