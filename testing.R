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
# when x22 and x23 don't include all x1's, gammas get weird
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

### get the start values once
start = hurdleIV.start_vals(regs, family = 'lognormal')
start = start[-5] # get rid of rho
start = tagBeg(start) # tag the beginning of the pi values
likelihood = getLik('lognormal')
out = optim(start,loglik_lgnorm)

### do it from true values
start = lgiv_params_mult[-c(1,2,3,4)]
start = tagBeg(start) # tag the beginning of the pi values
likelihood = getLik('lognormal')
out = optim(start,loglik_lgnorm)



