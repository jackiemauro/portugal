###############################################################
##
## A function which outputs start values for a hurdle IV 
## The input should be of the form:
#
# regs = list(formula = y1 ~ x11 + x21 + x22 + x23,
#             exog = list(x11,x12),
#             inst = list(z1,z2,z3),
#             endog = list(x21,x22,x23),
#             start_val = F,
#             endogReg = list(x21~x11+x12+z1+z2+z3,
#                          x22~x11+z1+z2+z3,
#                          x23~x11+z1+z2+z3))
# 
## set y1 to the log of y1 if you want to run a loglinear model
## if you leave endogReg empty, will default to using all exog & inst
## in the first stage regressions
###############################################################

# the start values function should run (generalized) linear regressions
# it should return two for y1 and one for each endogenous regressor
hurdleIV.start_vals <- function(object,family,...){
  
  if(family == 'lognormal'){
    formula = update(as.formula(object$formula), log(y1)~.)
  }
  else{
    formula = as.formula(object$formula)
  }
  
  # get regression for y1>0
  linReg = lm(formula, subset = y1>0)
  
  # get probit regression
  temp = as.formula(object$formula)
  formula = update.formula(temp, as.formula(paste('I(',temp[[2]],'>0)~.')))
  probReg = glm(formula, family = binomial(link = "probit"))
  
  # get endogenous regressions
  # if the user has specified the regressions, use those
  if(length(object$endogReg)>0){
    endogList = object$endogReg
    endRegs = list()
    for(i in 1:length(endogList)){
      formula = endogList[[i]]
      res = lm(formula)
      endRegs[[i]] = res
    }
  }
  # if the user has not specified, put all exogenous variables in 
  # all first stage regressions
  else{
    exo = matrix(unlist(object$exog),ncol = length(object$exog),byrow=T)
    instr = matrix(unlist(object$inst),ncol = length(object$inst),byrow = T)
    x = cbind(exo,instr)
    endRegs = list()
    for(i in 1:length(object$endog)){
      y = object$endog[[i]]
      endRegs[[i]] <- lm(y~x)
    }
  }
  
  j = length(object$endog); l = length(object$inst); m = length(linReg$coefficients)
  beta1 = linReg$coefficients[1:(m-j)]; beta2 = linReg$coefficients[(m-j+1):m]
  sig_u = sd(linReg$res)
  gamma1 = probReg$coefficients[1:(m-j)]; gamma2 = probReg$coefficients[(m-j+1):m]
  pi2 = list(); pi1 = list(); sig_v = list()
  for(i in 1:j){
    m = length(endRegs[[i]]$coefficients)
    pi1[[i]] = endRegs[[i]]$coefficients[1:(m-l)]
    pi2[[i]] = endRegs[[i]]$coefficients[(m-l+1):m]
    sig_v[[i]] = sd(endRegs[[i]]$res)
  }
  
  out = list(sig_u = sig_u, sig_v = sig_v,
             tau0 = rep(0,j), tau1 = rep(0,j), rho = 0,
             beta1 = beta1, beta2 = beta2,
             gamma1 = gamma1, gamma2 = gamma2,
             pi1 = pi1, pi2 = pi2)
  out
}
