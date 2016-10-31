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
source("loglik_max.R")
source("randStartFn.R")

require(MASS)
require(truncnorm)
require(xtable)

### get the start values and run
full.run <- function(regs
                     , method = 'Nelder-Mead'
                     , maxit = 5000
                     , loglik = loglik_lgnorm
                     , trace = 1){
  
  start = hurdleIV.start_vals(regs, family = 'lognormal')
  start = start[-5] # get rid of rho
  start = tagBeg(start) # tag the beginning of the pi values
  #likelihood = getLik('lognormal')
  a = proc.time()
  out = optim(start
              , loglik
              , method= method
              , hessian  = T
              , control = list(maxit = maxit,trace = trace)
  )
  time = proc.time() - a
  
  res = data.frame(true = tagBeg(pars)
                   , est = out$par
                   , start = start
                   , vars = diag(out$hessian))
  
  lls = data.frame(True = loglik_lgnorm(tagBeg(pars))
                   , Est = loglik_lgnorm(out$par)
                   , Reg = loglik_lgnorm(start))
  
  output = list(Results = res
                , Logliks = lls
                , Hessian = out$hessian
                , Gradient = out$counts['gradient']
                , time = time)
  
  print(paste('max eigenvalue = ', max(eigen(out$hessian)$value)))
  print(paste('min eigenvalue = ', min(eigen(out$hessian)$value)))
  return(output)
}


# output histograms
see.dist <- function(data=toprint,x,sub = F,eigs = F){
#   if(sub == T){data = data[which(data[,15]==0),]}
#   if(eigs == T){data = data[which(data[])]}
  m = dim(data)[1]-1
  hist(as.numeric(data[1:m,x])
       ,main = nm[x]
       ,xlab = 'n')
  abline(v = as.numeric(data[m+1,x]), col = 'red')
  abline(v = mean(as.numeric(data[1:m,x]),na.rm = T), col = 'black',lwd = 2)
  legend('topright'
         , col = c('red', 'black')
         , c('truth', 'average')
         , lty = c(1,1)
         , lwd = c(1,2)
  )
}


########### 1-1-1 formation ##########
regs = list(formula = 'y1 ~ x11 + x21',
            exog = list('x11'),
            inst = list('z1'),
            endog = list('x21'),
            start_val = F,
            endogReg = list('x21~x11+z1'))

detach(data)
sims = hurdleIV.gen_hurdleSim(n = 5000,
                              family = 'lognormal'
                              ,params = lgiv_params
                              ,formula = regs)

data = sims$dat
pars = sims$parameters
attach(data)

plot(x11~x21)

output = full.run(regs)

# compare to max function -- should be the same
mpar = c(pars$sig_v, pars$sig_u, pars$tau1, pars$tau0
         , pars$beta1, pars$beta2, pars$gamma1, pars$gamma2
         , pars$pi1[[1]][1], pars$pi2[[1]][1])
x2 = x21; x1 = x11; z = z1
loglik_max(mpar)
loglik_lgnorm(tagBeg(pars))

# # do random starts
# randstarts111 = rand.Start(100)
# 
# # output results
# write.csv(randstarts111$startvals, file = "RandStart111.csv")
# 
# trll = loglik_lgnorm(tagBeg(pars))
# truth = c(unlist(pars),trll,NA,NA)
# toprint = rbind(cbind(randstarts111$estimates
#                       ,randstarts111$loglikelihood
#                       ,randstarts111$convergence
#                       ,apply(randstarts111$eigenvalues, 1, min))
#                 ,truth)
# nm111 = c(names(pars), 'll', 'conv', 'eig')
# colnames(toprint)<-nm111
# write.csv(toprint, "randOut111.csv")

rs111 = read.csv('randOut111.csv')
rs111 = rs111[,2:dim(rs111)[2]]
names(rs111)[(dim(rs111)[2]-2):dim(rs111)[2]] <- c('ll', 'conv', 'mineign')
conv = which(rs111$conv == 0)
goodEig = which(rs111$mineign > 0)
conv_good = which(rs111$conv == 0 & rs111$mineign > 0)

for(i in 1:(dim(rs111)[2] - 2)){
  see.dist(data = rs111, x = i)
}

# look only at converged
for(i in 1:(dim(rs111)[2] - 2)){
  see.dist(data = rs111[conv,], x = i)
}


# look only at good eig's
for(i in 1:(dim(rs111)[2] - 2)){
  see.dist(data = rs111[goodEig,], x = i)
}

# look only at the really good ones
for(i in 1:(dim(rs111)[2] - 2)){
  see.dist(data = rs111[conv_good,], x = i)
}

hist(rs111$ll[conv_good])


########## 1-1-2 formation ##########
detach(data)

regs = list(formula = 'y1 ~ x11 + x12 + x21',
            exog = list('x11', 'x12'),
            inst = list('z1'),
            endog = list('x21'),
            start_val = F,
            endogReg = list('x21~x11+x12+z1'))

# get rid of coefficients on variables we don't have anymore
pars_new = lgiv_params_mult
pars_new$beta2 = pars_new$beta2[1]
pars_new$gamma2 = pars_new$gamma2[1]
pars_new$sig_v = pars_new$sig_v[1]
pars_new$tau1 = pars_new$tau1[1];pars_new$tau0 = pars_new$tau0[1]
pars_new$pi1 = list(pars_new$pi1[[1]])
pars_new$pi2 = list(pars_new$pi2[[1]][1])
pars_new$muz = pars_new$muz[1]
pars_new$sigz = 2

sims = hurdleIV.gen_hurdleSim(n = 5000,
                              family = 'lognormal', 
                              params = pars_new,
                              formula = regs)

data = sims$dat
pars = sims$parameters
attach(data)

par(mfrow = c(1,3))
plot(x11~x12); plot(x11~x21); plot(x12~x21)
par(mfrow = c(1,1))

output = full.run(regs)
# bad hessian if sigz too low (0.1)
# bump it up to 2 and all eig's positive
# still crazy high variance

View(output$Results)

betaH = result$hessian[5:7, 5:7]; eigen(betaH)$va
gammaH = result$hessian[8:10, 8:10]; eigen(gammaH)$va
piH = result$hessian[11:13, 11:13]; eigen(piH)$va

# # run at random start points
# randStart112 = rand.Start(100)
# 
# # output results
# write.csv(randStart112$startvals, file = "RandStart112.csv")
# 
# 
# trll = loglik_lgnorm(tagBeg(pars))
# truth = c(unlist(pars),trll,NA,NA)
# toprint = rbind(cbind(randStart112$estimates
#                       ,randStart112$loglikelihood
#                       ,randStart112$convergence
#                       ,apply(randStart112$eigenvalues, 1, min))
#                 ,truth)
# nm<-c(names(start),'likelihod','unconverged','min eigenval')
# write.csv(toprint, "randOut112.csv")

rs112 = read.csv('randOut112.csv')
rs112 = rs112[,2:dim(rs112)[2]]
names(rs112)[(dim(rs112)[2]-2):dim(rs112)[2]] <- c('ll', 'conv', 'mineign')
conv = which(rs112$conv == 0)
goodEig = which(rs112$mineign > 0)

for(i in 1:(dim(rs112)[2] - 2)){
  see.dist(data = rs112, x = i)
}

# look only at converged
for(i in 1:(dim(rs112)[2] - 2)){
  see.dist(data = rs112[conv,], x = i)
}


# look only at good eig's
for(i in 1:(dim(rs112)[2] - 2)){
  see.dist(data = rs112[goodEig,], x = i)
}


########### n-n-n formation ##############
regs = list(formula = 'y1 ~ x11 + x12+ x21 + x22 + x23',
            exog = list('x11','x12'),
            inst = list('z1','z2','z3'),
            endog = list('x21','x22','x23'),
            start_val = F,
            endogReg = list('x21~x11+x12+z1+z2+z3',
                            'x22~x11+z1+z2+z3',
                            'x23~x12+z1+z2+z3'))
detach(data)
sims = hurdleIV.gen_hurdleSim(n = 5000,
                              family = 'lognormal'
                              ,params = lgiv_params_mult
                              #,params = lgiv_params
                              ,formula = regs)

data = sims$dat
pars = sims$parameters
attach(data)

### do it from true values
# pretty good
start = pars
start = tagBeg(start) # tag the beginning of the pi values
likelihood = getLik('lognormal')
outTrue = optim(start,loglik_lgnorm)
View(cbind(tagBeg(start),outTrue$par))
llTrue = loglik_lgnorm(start)
llEst = loglik_lgnorm(outTrue$par)



### do it from true values
# pretty good
start = pars
start = tagBeg(start) # tag the beginning of the pi values
likelihood = getLik('lognormal')
outTrue = optim(start,loglik_lgnorm)
View(cbind(tagBeg(start),outTrue$par))
llTrue = loglik_lgnorm(start)
llEst = loglik_lgnorm(outTrue$par)

### get the start values and run
# not great
start = hurdleIV.start_vals(regs, family = 'lognormal')
start = start[-5] # get rid of rho
start = tagBeg(start) # tag the beginning of the pi values
likelihood = getLik('lognormal')
a = proc.time()
out = optim(start
            , likelihood
            , method= "Nelder-Mead"
            , hessian  = T
            , control = list(maxit = 7000,trace = 1)
            )
timeSoLongUgh = proc.time() - a

View(cbind(tagBeg(pars),out$par,start,diag(out$hessian)))
llTrue = loglik_lgnorm(tagBeg(pars))
llEst = loglik_lgnorm(out$par)
llReg = loglik_lgnorm(start)

max(eigen(out$hessian)$value); min(eigen(out$hessian)$value)
out$counts


# Max q's
# a) If you do lots of random starting points, do you always converge to the same bad point?
n = 100
randOut = NA; llEst = NA; llStart = NA
randStart = matrix(c(rep(NA,n*length(start))),ncol = length(start))

for(i in 1:n){
  
  print(paste("round ",i," of ",n, sep = ""))
  val = sapply(start, function(x) x + rnorm(1,0,sd = abs(x*2)))
  result = tryCatch({
    optim(val,loglik_lgnorm,control = list(maxit=5000))
    },error = function(e){
      print(paste("optimizer error, skipping round ",i,sep = ""))
      return(NA)
    }
    )
#   if(!is.na(result[1])){
#     llEst[i] = loglik_lgnorm(result)
#     llStart[i] = loglik_lgnorm(val)
#   }
#   else{
#     llEst[i] = NA
#     llStart[i] = NA
#   }
  randStart[i,] = val
  randOut[i] <- result
}
LS.df = as.data.frame(do.call(cbind, randOut))
test = LS.df[, colSums(is.na(LS.df)) != nrow(LS.df)]
lls = apply(test,2,loglik_lgnorm)
test = rbind(test,lls)
rownames(test)[34]<-"log likelihood"
write.csv(test, "randStartoutput_moreVar.csv")


write.csv(t(randStart),"randStartPoints_moreVar.csv")

# b) What does the hessian look like at the point you converge to?  
# Is it well-conditioned?
start = hurdleIV.start_vals(regs, family = 'lognormal')
start = start[-5] # get rid of rho
start = tagBeg(start) # tag the beginning of the pi values
likelihood = getLik('lognormal')
out = optim(start,likelihood, hessian  = T)

llTrue = loglik_lgnorm(tagBeg(pars))
llEst = loglik_lgnorm(out$par)
llReg = loglik_lgnorm(start)

View(cbind(c(llTrue,round(tagBeg(pars),4)),
           c(llEst,round(out$par,4)),
           c(llReg,round(start,4))))

hessRatio = max(eigen(out$hessian)$values)/min(eigen(out$hessian)$values)


# c) What is the value of the objective at these minima versus the true value?  Are they much lower?
#see above


######### checking gradient (BFGS) ############

out = optim(start
            , likelihood
            , method= "BFGS"
            , hessian  = T
            , control = list(maxit = 5000,trace = 1)
)
max(out$hessian); min(out$hessian)
out$counts['gradient']

out = optim(start
            , likelihood
            , method= "Nelder-Mead"
            , hessian  = T
            , control = list(maxit = 5000,trace = 1)
)



########## more instruments than endogenous variables ##########
# does not work from estimated start values
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