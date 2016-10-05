### no instrument

beta1 = NA; beta2 = NA; sd = NA
for(i in 1:1000){
  x = rnorm(1000,3,4)
  z = rnorm(1000,5,2)
  
  Sig = mvrnorm(1000, c(0,0,0), 
                 matrix(c(2,.3,-.6,
                          .3,5,.4,
                          -.6,.4,2), 
                        byrow = T, ncol = 3))
  
  x2 = 4*z + x + Sig[,3]
  y1 = 5*x + 3*x2 +  Sig[,1]
  y2 = -3*x + 2*x2 + Sig[,2]
  y = y1 + y2
  
  r = lm(y ~ x + x2)
  beta1[i] = r$coef[2]; beta2[i] = r$coef[3]
  sd[i] = sd(r$resid) 
}


### simple iv reg
beta2 = NA
for(i in 1:1000){
  x = rnorm(1000,3,4)
  z = rnorm(1000,5,2)
  
  Sig = mvrnorm(1000, c(0,0,0), 
                matrix(c(2,.3,-.6,
                         .3,5,.4,
                         -.6,.4,2), 
                       byrow = T, ncol = 3))
  
  x2 = 4*z + x + Sig[,3]
  y1 = 5*x + 3*x2 +  Sig[,1]
  y2 = -3*x + 2*x2 + Sig[,2]
  y = y1 + y2
  
  r1 = lm(x2 ~ x + z)
  pred1 = predict(r1)
  full1 = lm(y1 ~ x + pred1)
  full2 = lm(y2 ~ x + pred1)
  full = lm(y ~ x + pred1)
  
  beta2[i] = full$coef[3]
}

### with full fn
source("startValuesFn.R")
source("loglik_lgnorm.R")
source("parameter_values.R")
source("hurdleIVClass.R")
source("dataSimFn.R")

require(MASS)
require(truncnorm)
require(xtable)

regs = list(formula = 'y1 ~ x11 + x21 ',
            exog = list('x11'),
            inst = list('z1'),
            endog = list('x21'),
            start_val = F,
            endogReg = list('x21~x11+z1'))

pars = lgiv_params
pars$pi1 = list(c(.1,2)); pars$pi2 = list(.7)

detach(data)
beta2 = NA; beta21 = NA; beta22 = NA
for(i in 1:10){
  
  #simulate dataset
  sims = hurdleIV.gen_hurdleSim(n = 1000,
                                family = 'lognormal', 
                                params = pars,
                                formula = regs)
  
  data = sims$dat
  params = sims$parameters
  
  #decompose y1
  ## one option:
  y11 = exp(log(data$y1) * .3 ); y11[data$y1==0] = 0
  y12 = exp( log(data$y1) - log(y11) ); y12[data$y1==0] = 0
  
  ## more complicated fn:
  #y11 = exp(log(data$y1)*.4 - .1*(log(data$y1))^2 ); y11[data$y1==0] = 0
  #y12 = exp( log(data$y1) - log(y11) ); y12[data$y1==0] = 0

  # run the regressions
  attach(data)
  start = params
  start = tagBeg(start) # tag the beginning of the pi values
  likelihood = getLik('lognormal')
  out = optim(start,likelihood)
  beta2[i] = out$par['beta2_ls6.elem1']
  
  detach(data)
  data$y1 = y11
  attach(data)
  out1 = optim(start,likelihood)
  beta21[i] = out1$par['beta2_ls6.elem1']
  
  detach(data)
  data$y1 = y12
  attach(data)
  out2 = optim(start,likelihood)
  beta22[i] = out2$par['beta2_ls6.elem1']
  detach(data)
}


mean(beta21+beta22-beta2)