# is probit then lognormal consistent?
# looks like it is. 
# our approach does marginally better estimating coefficients, but not much.
rm(list = ls())
source("wrapperFunc.R")
source("hurdleIVsim.R")


# choose your gamma and beta values
gammas = seq(-5,5,length=10)
beta2 = .02
rhos = seq(0,.5,length=10)

# make your datasets
make.dat.gamma <- function(gamma2){
  hurdle.IV.sim(gamma = c(.2,.4,gamma2)
                ,beta = c(.05,-.1,beta2)
                ,type = "cragg",silent = T)
}

make.dat.rho <- function(rho){
  hurdle.IV.sim(rho = rho
                ,beta = c(.05,-.1,beta2)
                ,type = "cragg",silent = T)
}


# naive cragg function (runs ivprobit then ivreg)
naive_cragg <- function(dat){
  require('AER')
  require('ivprobit')
  lin = lm(endog~exog1 + inst1, data = dat)
  preds = predict(lin)
  prob = ivprob(y=dat$y0,x1=dat$exog1,y2=dat$endog,x=data.frame(dat$exog1,dat$inst1))
  reg = ivreg(y ~ exog1 + endog | exog1 + inst1,data = dat, subset = y>0)
  beta = reg$coef['endog']
  gamma = prob$coeff[3]
  return(beta)
}

# test for changing values of gamma
n = 100
out = matrix(rep(NA,n*length(gammas)),ncol = length(gammas))
for(i in 1:n){
  out[i,] = unlist(lapply(as.list(gammas), function(x) naive_cragg(make.dat.gamma(x))))
}

means = apply(out,2,mean)
sds = apply(out,2,sd)

#png(filename = "naive_cragg.png")
plot(gammas,means, xlab = expression(gamma), ylab = expression(beta)
     ,ylim = c(min(means-2*sds),max(means+2*sds)))
abline(beta2,0)
arrows(gammas, means-2*sds, gammas, means+2*sds, length=0.05, angle=90, code=3)
text(beta2*2,expression("True " * beta))
#dev.off()

# test for changing values of rho -- something is wrong
n = 10
b = matrix(rep(NA,n*length(rhos)),ncol = length(rhos))
for(i in 1:n){
  b[i,] = unlist(lapply(as.list(rhos), function(x) naive_cragg(make.dat.rho(x))))
}

means2 = apply(b,2,mean)
sds2 = apply(b,2,sd)

png(filename = "naive_cragg_chg_rho.png")
plot(rhos,means2, xlab = expression(rho), ylab = expression(beta)
     ,ylim = c(min(means2-2*sds2),max(means2+2*sds2)))
abline(beta2,0)
arrows(rhos, means2-2*sds2, rhos, means2+2*sds2, length=0.05, angle=90, code=3)
text(beta2*2,expression("True " * beta))
dev.off()

