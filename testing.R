#####
# trash/testing file
####

setwd("C:/Users/jackie/Dropbox/portugal/theory paper")

rm(list = ls())
par(cex = 1.25)
par.old <- par()

source("startValuesFn.R")
source("loglik_lgnorm.R")
source("loglik_cragg.R")
source("parameter_values.R")
source("hurdleIVClass.R")
source("dataSimFn.R")

require(MASS)
require(truncnorm)
require(xtable)

# reg = hurdleIV(y1 ~ x11 + x21 + x22 + x23,
#              exog = list(x11,x12),
#              endog = list(x21~x11+x12+z1+z2+z3,
#                           x22~x11+z1+z2+z3,
#                           x23~x12+z1+z2+z3))


# for now, need to set all these values to strings
regs = list(formula = 'y1 ~ x11 + x12+ x21 + x22 + x23',
            exog = list('x11','x12'),
            inst = list('z1','z2','z3'),
            endog = list('x21','x22','x23'),
            start_val = F,
            endogReg = list('x21~x11+x12+z1+z2+z3',
                            'x22~x11+z1+z2+z3',
                            'x23~x12+z1+z2+z3'))

detach(data)
sims = hurdleIV.gen_hurdleSim(n = 1000,
                              family = 'lognormal', 
                              params = lgiv_params_mult,
                              formula = regs)

data = sims$dat
pars = sims$parameters
attach(data)

########### running the code ##############

### get the start values and run
# not great
start = hurdleIV.start_vals(regs, family = 'lognormal')
start = start[-5] # get rid of rho
start = tagBeg(start) # tag the beginning of the pi values
likelihood = getLik('lognormal')
out = optim(start,likelihood)

### do it from true values
# pretty good
start = pars
start = tagBeg(start) # tag the beginning of the pi values
likelihood = getLik('lognormal')
out = optim(start,likelihood)

### do it with cragg -- not really functional except from true
start = hurdleIV.start_vals(regs, family = 'cragg1')
start = start[-5] # get rid of rho
start = tagBeg(start) # tag the beginning of the pi values
likelihood = getLik('cragg1')
out = optim(start,likelihood)

### from true values
start = cr_params_mult[-c(1,2,3,4)]
start = tagBeg(start) # tag the beginning of the pi values
likelihood = getLik('cragg1')
out = optim(start,likelihood)

### more instruments than endogenous variables
# does not work
detach(data)
regs = list(formula = 'y1 ~ x11 + x21',
            exog = list('x11'),
            inst = list('z1','z2'),
            endog = list('x21'),
            start_val = F,
            endogReg = list('x21~x11+z1+z2'))

# get rid of coefficients on variables we don't have anymore
pars_new = lgiv_params_mult
pars_new$beta2 = pars_new$beta2[1]
pars_new$gamma2 = pars_new$gamma2[1];pars_new$sig_v = pars_new$sig_v[1]
pars_new$tau1 = pars_new$tau1[1];pars_new$tau0 = pars_new$tau0[1]
pars_new$pi1 = list(pars_new$pi1[[1]][c(1,2)])
pars_new$pi2 = list(pars_new$pi2[[1]][c(1,2)])
#pars_new$mux1 = pars_new$mux1[1];
pars_new$mux1 = 1.3
#pars_new$sigx1 = pars_new$sigx1[1]
pars_new$sigx1 = c(0.1,0.1)
pars_new$muz = pars_new$muz[-3]; 
#pars_new$sigz = pars_new$sigz[-3]
pars_new$sigz = c(0.1,0.1)
pars_new$beta1 = pars_new$beta1[-3]; pars_new$gamma1 = pars_new$gamma1[-3]


sims = hurdleIV.gen_hurdleSim(n = 100,
                              family = 'lognormal', 
                              params = pars_new,
                              formula = regs)

detach(data)
data = sims$dat
pars = sims$parameters
attach(data)

### get the start values and run
#bad
start = hurdleIV.start_vals(regs, family = 'lognormal')
start = start[-5] # get rid of rho
start = tagBeg(start) # tag the beginning of the pi values
likelihood = getLik('lognormal')
out = optim(start,loglik_lgnorm)
View(cbind(out$par,unlist(pars)))

### do it from true
start = pars
start = tagBeg(start) # tag the beginning of the pi values
likelihood = getLik('lognormal')
out = optim(start,likelihood)
View(cbind(out$par,unlist(pars)))