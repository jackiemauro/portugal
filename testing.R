#####
# trash/testing file
####

setwd("C:/Users/jackie/Dropbox/portugal/theory paper")

rm(list = ls())
par(cex = 1.25)
par.old <- par()

source("startValuesFn.R")
source("loglik_lgnorm.R")
source("JM methods paper/loglik_cragg.R")
source("JM methods paper/parameter_values.R")
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

sims = hurdleIV.gen_hurdleSim(n = 1000,
                              family = 'lognormal', 
                              params = lgiv_params_mult,
                              formula = regs)

data = sims$dat
attach(data)

j = length(regs$endogReg); l = length(regs$inst)

########### running the code ##############

### get the start values and run
start = hurdleIV.start_vals(regs, family = 'lognormal')
start = start[-5] # get rid of rho
start = tagBeg(start) # tag the beginning of the pi values
likelihood = getLik('lognormal')
out = optim(start,likelihood)

### do it from true values
start = lgiv_params_mult[-c(1,2,3,4)]
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
detach(data)
regs = list(formula = 'y1 ~ x11 + x12+ x21 + x22',
            exog = list('x11','x12'),
            inst = list('z1','z2','z3'),
            endog = list('x21','x22'),
            start_val = F,
            endogReg = list('x21~x11+x12+z1+z2+z3',
                            'x22~x11+z1+z2+z3'))

# get rid of one set of coefficients
pars_new = lgiv_params_mult
pars_new$beta2 = pars_new$beta2[-3]
pars_new$gamma2 = pars_new$gamma2[-3];pars_new$sig_v = pars_new$sig_v[-3]
pars_new$tau1 = pars_new$tau1[-3];pars_new$tau0 = pars_new$tau0[-3]
pars_new$pi1 = pars_new$pi1[-3]; pars_new$pi2 = pars_new$pi2[-3]


sims = hurdleIV.gen_hurdleSim(n = 100,
                              family = 'lognormal', 
                              params = pars_new,
                              formula = regs)

data = sims$dat
attach(data)

### get the start values and run
start = hurdleIV.start_vals(regs, family = 'lognormal')
start = start[-5] # get rid of rho
start = tagBeg(start) # tag the beginning of the pi values
likelihood = getLik('lognormal')
out = optim(start,likelihood)