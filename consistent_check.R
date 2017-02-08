# is probit then lognormal consistent?
# looks like it is. 
# our approach does marginally better estimating coefficients, but not much.
rm(list = ls())
source("wrapperFunc.R")
source("hurdleIVsim.R")
source("hurdleIVsim_rawErrors.R")
source("startFn.R")
source("loglikFunc.R")
require('AER')
require('ivprobit')

beta2 = .02; gamma2 = -.07
beta = list(); gamma = list(); pi = list()
beta_ll = list(); gamma_ll = list(); pi_ll = list(); cov_ll = list()
for(i in 1:100){
  dat = hurdle.IV.sim_rawError(pi = c(1,-1,3), #intercept, exog, inst
                      gamma = c(.2,.4,gamma2), #intercept, exog, endog
                      beta = c(.05,-.1,beta2)  #intercept, exog, endog
                      )
  
  lin = lm(endog~exog1 + inst1, data = dat)
  pi[[i]] = lin$coef
  preds = predict(lin)
  prob = ivprob(y=dat$y0,x1=dat$exog1,y2=dat$endog,x=data.frame(dat$exog1,dat$inst1))
  gamma[[i]] = prob$coef
  lognorm = ivreg(log(y) ~ exog1 + endog | exog1 + inst1,data = dat, subset = y>0)
  beta[[i]] = lognorm$coef
  
  out = hurdle.IV(formula = y~endog + exog1
                  ,inst = inst1
                  ,endog = endog
                  ,exog = exog1
                  ,data = dat
                  ,endog_reg = list()
                  ,start_val = list()
                  )
  beta_ll[[i]] = out[[2]]$beta
  gamma_ll[[i]] = out[[2]]$gamma
  pi_ll[[i]] = out[[2]]$pi[[1]]
  cov_ll[[i]] = out[[2]]$cov
}

betas = matrix(unlist(beta), ncol = 3,byrow = T)
gammas = matrix(unlist(gamma), ncol = 3, byrow = T)
pis = matrix(unlist(pi), ncol = 3, byrow = T)

betas_ll = matrix(unlist(beta_ll), ncol = 3,byrow = T)
gammas_ll = matrix(unlist(gamma_ll), ncol = 3, byrow = T)
pis_ll = matrix(unlist(pi_ll), ncol = 3, byrow = T)

par(mfrow = c(1,2))
hist(betas[,3], main = "Sep beta", xlab = paste("SD: ",round(sd(betas[,3]),3)))
abline(v = mean((betas[,3])), col = "green")
abline(v = beta2, col = "red")

hist(betas_ll[,2], main = "FIML beta", xlab = paste("SD: ",round(sd(betas_ll[,2]),3)))
abline(v = mean((betas_ll[,2])), col = "green")
abline(v = beta2, col = "red")
par(mfrow = c(1,1))


par(mfrow = c(1,2))
hist(gammas[,3], main = "Sep gamma", xlab = paste("SD: ",round(sd(gammas[,3]),3)))
abline(v = mean((gammas[,3])), col = "green")
abline(v = gamma2, col = "red")

hist(gammas_ll[,2], main = "FIML gamma", xlab = paste("SD: ",round(sd(gammas_ll[,2]),3)))
abline(v = mean((gammas_ll[,2])), col = "green")
abline(v = gamma2, col = "red")
par(mfrow = c(1,1))
