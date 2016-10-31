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
pars$pi1 = list(c(2)); pars$pi2 = list(.7)

detach(data)
beta2 = NA; beta21 = NA; beta22 = NA; diff = NA
for(i in 1:50){
  
  #simulate dataset -- with or without intercept
  sims = hurdleIV.gen_hurdleSim(n = 1000,
                                family = 'lognormal', 
                                params = pars,
                                formula = regs)
  
  data = sims$dat
  params = sims$parameters
  
  #decompose y1
  ## one option:
  #y11 = exp(log(data$y1) * .3 ); y11[data$y1==0] = 0
  #y12 = exp( log(data$y1) - log(y11) ); y12[data$y1==0] = 0
  
  ## more complicated fn:
  y11 = exp(log(data$y1)*.4 - .1*(log(data$y1))^2 ); y11[data$y1==0] = 0
  y12 = exp( log(data$y1) - log(y11) ); y12[data$y1==0] = 0

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
  
  diff[i] = beta21[i]+beta22[i]-beta2[i]
}

mean(beta21+beta22-beta2)
x = seq(1,50,1)
plot(diff~x)
abline(0,0)


######## full approach #####
x = rnorm(1000,1.5,2)
z = rnorm(1000,1.3,2)

corr = 0 #if this is too high, matrix not PSD
Sig_err = matrix(c(1,0,0, .02,
                   0,.75,corr,.04,
                   0,corr,.3,.07,
                   .02,.4,.07,.8),
                 ncol = 4, byrow = T)
beta11 = .3; beta21 = -.4
beta12 = .2; beta22 = -.1
gamma2 = .5


A = matrix(c(1,0,0,gamma2,
             0,1,0,beta21,
             0,0,1,beta22,
             0,0,0,1),
           ncol = 4, byrow = T)

Sig = A%*%Sig_err%*%t(A)

  
euv = mvrnorm(1000,c(0,0,0,0),Sig)
eta = euv[,1]
u1 = euv[,2]
u2 = euv[,3]
v = euv[,4]

x2 = x*.4 + z*.7 + v
lny11 = x*beta11 + x2*beta21 + u1
lny12 = x*beta12 + x2*beta22 + u2
lny1 = lny11 + lny12

y0star = x*-.05 + x2*gamma2 + eta
y0 = as.numeric(y0star>0)

y = exp(lny1)*y0
y11 = exp(lny11)*y0
y12 = exp(lny12)*y0

# with real regression
data1 = data.frame(x11=x, x21 = x2, z1 = z, y1 = y11)
data2 = data.frame(x11=x, x21 = x2, z1 = z, y1 = y12)
data =  data.frame(x11=x, x21 = x2, z1 = z, y1 = y)

attach(data)
start = list(sig_u = .75, sig_v =  .8, tau0 = .02, tau1 = .04,
             beta1 = .5, beta2 = -.5, gamma1 = -.05, gamma2 = .5,
             pi1 = list(.2), pi2 = list(.7))
start = tagBeg(start)
out = optim(start,loglik_lgnorm)
beta2_res = out$par['beta2_ls6.elem1']
detach(data)

attach(data1)
start = list(sig_u = .75, sig_v =  .8, tau0 = .02, tau1 = .04,
             beta1 = beta11, beta2 = beta21, gamma1 = -.05, gamma2 = .5,
             pi1 = list(.2), pi2 = list(.7))
start = tagBeg(start)
out = optim(start,loglik_lgnorm)
beta21_res = out$par['beta2_ls6.elem1']
detach(data1)

attach(data2)
start = list(sig_u = .75, sig_v =  .8, tau0 = .02, tau1 = .04,
             beta1 = beta12, beta2 = beta22, gamma1 = -.05, gamma2 = .5,
             pi1 = list(.2), pi2 = list(.7))
start = tagBeg(start)
out = optim(start,loglik_lgnorm)
beta22_res = out$par['beta2_ls6.elem1']
detach(data2)

out = beta22_res + beta21_res - beta2_res 
