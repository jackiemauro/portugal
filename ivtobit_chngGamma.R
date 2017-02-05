### figure with iv tobit coefficients ##
rm(list = ls())
setwd("C:/Users/jackie/Desktop/own research/portugal/code attempt jan 2017")


# not giving the s-shape.
require(censReg)
source("wrapperFunc.R")
source("hurdleIVsim.R")
source("startFn.R")
source("loglikFunc.R")

parameters = list(pi = c(0,2,.7), #intercept, exog, inst
                  gamma = c(0,-.2,-.5), #intercept, exog, endog
                  beta = c(0,-.5,-.8), #intercept, exog, endog
                  endog_reg = list(),
                  exog_mean = .2,
                  exog_sd = sqrt(.2),
                  z_mean = 1.5,
                  z_sd = sqrt(.5),
                  endog_sd = 1.8,
                  y_sd = 1.2,
                  rho = 0,
                  tau0 = .04,
                  tau1 = .02)

#tobit assumption values
newpar = parameters
newpar$tau1 = 0
newpar$rho = 0
newpar$y_sd <- 1
sig_y0_x2 <- 1 - parameters$tau0^2/parameters$endog_sd^2
assmp_b2 <- (parameters$gamma[3] - parameters$tau0/parameters$endog_sd^2)/sig_y0_x2
assmp_b1 <- parameters$gamma[1:2]/sig_y0_x2
assmp_g1 <- parameters$beta[1:2]*sig_y0_x2
assmp_g2 <- parameters$beta[3]*sig_y0_x2 + parameters$tau0/parameters$endog_sd^2
newpar$gamma = c(assmp_g1, assmp_g2)

dat = hurdle.IV.sim(formula = F,
                    n = 1000,
                    pi = newpar$pi,
                    gamma = newpar$gamma,
                    #gamma = c(newpar$gamma[1:2],-1),
                    beta = newpar$beta,
                    exog_mean = newpar$exog_mean,
                    z_mean = newpar$z_mean,
                    exog_sd = newpar$exog_sd,
                    z_sd = newpar$z_sd,
                    endog_sd = newpar$endog_sd,
                    y_sd = newpar$y_sd,
                    rho = newpar$rho,
                    tau0 = newpar$tau0,
                    tau1 = newpar$tau1,
                    type = "cragg",
                    cond = F)

source("hurdleIVsim_badErrors.R")
dat = hurdle.IV.sim_badError(formula = F,
                    n = 1000,
                    pi = newpar$pi,
                    gamma = newpar$gamma,
                    #gamma = c(newpar$gamma[1:2],-1),
                    beta = newpar$beta,
                    exog_mean = newpar$exog_mean,
                    z_mean = newpar$z_mean,
                    exog_sd = newpar$exog_sd,
                    z_sd = newpar$z_sd,
                    endog_sd = newpar$endog_sd,
                    y_sd = newpar$y_sd,
                    rho = newpar$rho,
                    tau0 = newpar$tau0,
                    tau1 = newpar$tau1,
                    type = "cragg")

require(censReg)
reduced.form <- lm(endog~exog1 + inst1, dat)
attach(dat)
consistent.tobit <- censReg(y~fitted(reduced.form)+residuals(reduced.form) + exog1)
detach(dat)
coef(consistent.tobit)


gRange2 = seq(from = assmp_g2-2, to = assmp_g2+2, by = .25)

tobit_beta2 = rep(NA,length(gRange2))
sd = rep(NA,length(gRange2))
n = 1000; m = 100
for(ii in 1:length(gRange2)){
  
  print(paste("round ",ii))
  gamma = parameters$gamma
  gamma[3] = gRange2[ii]
  vals = c(rep(NA,m))
  for(jj in 1:m){
    # simulate data
    dat = hurdle.IV.sim(formula = F,
                        n = n,
                        pi = newpar$pi,
                        gamma = gamma,
                        beta = newpar$beta,
                        exog_mean = newpar$exog_mean,
                        z_mean = newpar$z_mean,
                        exog_sd = newpar$exog_sd,
                        z_sd = newpar$z_sd,
                        endog_sd = newpar$endog_sd,
                        y_sd = newpar$y_sd,
                        rho = newpar$rho,
                        tau0 = newpar$tau0,
                        tau1 = newpar$tau1,
                        type = "cragg",
                        silent = T)
    
    # iv tobit
    require(censReg)
    reduced.form <- lm(endog~exog1 + inst1, dat)
    attach(dat)
    consistent.tobit <- censReg(y~fitted(reduced.form)+residuals(reduced.form) + exog1)
    detach(dat)
    vals[jj] = coef(consistent.tobit)['fitted(reduced.form)']
  }
  tobit_beta2[ii] = mean(vals, na.rm = T)
  sd[ii] = sd(vals, na.rm = T)
}

dat_tob = data.frame(Gamma = gRange2, Beta = tobit_beta2, err = sd)
g <- ggplot(dat_tob,aes(x = Gamma, y = Beta)) +
  geom_point(size = 2)+
  geom_errorbar(aes(ymin=Beta-2*err, ymax=Beta+2*err), width=.1) +
  geom_hline(yintercept=assmp_b2, lty = 2)+
  geom_vline(xintercept=assmp_g2) +
  theme_bw()
g