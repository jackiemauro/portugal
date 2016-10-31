# for now, need to set all these values to strings
regs = list(formula = 'y1 ~ x11 + x21 ',
            exog = list('x11'),
            inst = list('z1'),
            endog = list('x21'),
            start_val = F,
            endogReg = list('x21~x11+z1'))

pars1 = lgiv_params
pars1$beta2 = .3; pars1$beta1[2] = -.5
pars1$pi1 = list(c(.1,2)); pars1$pi2 = list(.7)
sims1 = hurdleIV.gen_hurdleSim(n = 1000,
                               family = 'lognormal', 
                               params = pars1,
                               formula = regs)

data1 = sims1$dat
params1 = sims1$parameters

pars2 = pars1
pars2$beta2 = .5; pars2$beta1[2] = -.3
sims2 = hurdleIV.gen_hurdleSim(n = 1000,
                               family = 'lognormal', 
                               params = pars2,
                               formula = regs)

data2 = sims2$dat
params2 = sims2$parameters

y1 = data1$y1; y2 = data2$y1

lny = log(y1) + log(y2)


####### approach2

# works
r1 = lm(x2~ x1 + z)
pred1 = predict(r1)

lm1 = lm(lny1 ~ x1 + pred1)
lm2 = lm(lny2 ~ x1 + pred1)
lm = lm(lny ~ x1 + pred1)

lm$coef[3] - lm1$coef[3] - lm2$coef[3]

## use real fn and loop
Sig = matrix(c(.75,0,.04,0,1,.02,.04,.02,.8), byrow = T, ncol = 3)
beta11 = c(-.3); beta12 = c(.15)
beta21 = -.4; beta22 = -.1
params = lgiv_params
params$pi1 = list(c(.2)); params$pi2 = list(c(.7))

beta2_res = NA; beta21_res = NA; beta22_res = NA
for(i in 1:10){
  x1 = rnorm(1000,2,.6)
  z = rnorm(1000,.5,.5)
  
  euv = mvrnorm(1000,c(0,0,0),Sig)
  # think that this untransformed error term is a problem
  # then the question becomes -- which beta to use in the transformation?
  # and it seems odd that this should lead to bias in beta2
  v = euv[,3]; eta = euv[,2]; u = euv[,1]
  
  x2 = .1 + 2*x1 + .7*z + v
  y0star = .3 - .2*x1 + -.05*z + eta
  y0 = as.numeric(y0star>0)
  
  lny1 = beta11*x1 + beta21*x2 + .5*u
  lny2 = beta12*x1 + beta22*x2 + .5*u
  
  #lny = lny1 + lny2
  #lny = beta11[1] + beta12[1] + (beta11[2]+beta12[2])*x1 + (beta21 + beta22)*x2 + u
  lny = 0.15*x1 - .5*x2 + u
  
  y11 = exp(lny1); y11[y0==0] = 0
  y12 = exp(lny2); y12[y0==0] = 0
  y = exp(lny); y[y0==0] = 0
    
  # with real regression
  data1 = data.frame(x11=x1, x21 = x2, z1 = z, y1 = y11)
  data2 = data.frame(x11=x1, x21 = x2, z1 = z, y1 = y12)
  data =  data.frame(x11=x1, x21 = x2, z1 = z, y1 = y)
  
  attach(data)
  start = params
  start = tagBeg(start)
  out = optim(start,loglik_lgnorm)
  beta2_res[i] = out$par['beta2_ls8.elem1']
  detach(data)
  
  attach(data1)
  start = params
  start$beta2 = beta21
  start$beta1 = beta11
  start = tagBeg(start)
  out1 = optim(start,loglik_lgnorm)
  beta21_res[i] = out1$par['beta2_ls8.elem1']
  detach(data1)
  
  attach(data2)
  start = params
  start$beta2 = beta22
  start$beta1 = beta12
  start = tagBeg(start)
  out2 = optim(start,loglik_lgnorm)
  beta22_res[i] = out2$par['beta2_ls8.elem1']
  detach(data2)
}
