### calculates the log likelihood for a tobit IV regression
# does not include intercept terms

loglik_tobitiv <- function(t){
  
  # re-listify
  sig_u = t[grep('sig_u',names(t))]
  sig_v = t[grep('sig_v',names(t))]
  rho = t[grep('rho',names(t))]
  beta1 = t[grep('beta1',names(t))]
  beta2 = t[grep('beta2',names(t))]
  pi1 = t[grep('pi1',names(t))]
  pi2 = t[grep('pi2',names(t))]
  
  j = length(regs$endog)
  l = length(regs$inst)
  
  regStarts = c(grep('subelem1',names(pi1)),length(pi1)+1)
  v = diff(regStarts)
  pi1 = split(pi1, rep(1:length(v),v))
  
  regStarts = c(grep('subelem1',names(pi2)),length(pi2)+1)
  v = diff(regStarts)
  pi2 = split(pi2, rep(1:length(v),v))
  
  # make sig_err
  a = matrix(diag(j)*(unlist(sig_v))^2, nrow = j)
  b = cbind(rho,a)
  Sig_err = rbind(c(sig_u^2,rho),b)
  if(min(eigen(Sig_err)$values)<=0){
    print("negative eigenvalue in Sigma_err")
    return(Inf)
  }
  if((sig_u<=0)){
    print('negative variance value')
    return(Inf)
  }
  
  #Transformation matrix from (eta, u, v) to (y0*, log y1*, x2) (without the mean)
  A = diag(j+1)
  A[1,] = c(1,beta2)
  
  Sig = A%*%Sig_err%*%t(A)
  if(min(eigen(Sig)$values)<=0){
    print("negative eigenvalue in Sigma")
    return(Inf)
  }
  if(sig_u<=0){
    print('negative variance value')
    return(Inf)
  }
  
  # get the means
  censored = y1<=0
  
  # first get the x2 means so they're in the right shape
  mu_x2 = matrix(c(rep(NA,j*length(y1>0))),ncol = j)
  for(i in 1:j){
    formula = regs$endogReg[[i]]
    #mf = model.frame(formula = formula, na.action = NULL)
    mf = model.frame(formula = formula, na.action = NULL)
    m = dim(mf)[2]
    x <- model.matrix(attr(mf, "terms"), data=mf)
    x1temp = as.matrix(x[,2:(m-l)]); ztemp = as.matrix(x[,(m-l+1):m])
    mu_x2[,i] = x1temp%*%as.matrix(pi1[[i]]) + ztemp%*%as.matrix(pi2[[i]])
  }
  
  # then get the means for the y regressions
  formula = regs$formula
  #mf = model.frame(formula = formula, na.action = NULL)
  mf = model.frame(formula = formula)
  m = dim(mf)[2]
  x <- model.matrix(attr(mf,"terms"),data = mf)
  x1 = as.matrix(x[,2:(m-j)]); x2 = as.matrix(x[,(m-j+1):m])
  mu_y1 = x1%*%beta1 + mu_x2%*%beta2
  
  # now get the conditional means
  #Parameters for x2
  k = dim(Sig)[1]
  sig2_x2 = as.matrix(Sig[2,2:k])
  
  #Parameters for y1star given x2
  mu_y1_x2 = mu_y1 + c(t(Sig[1,2:k])%*%solve(Sig[2,2:k])%*%t(x2-mu_x2))
  sig2_y1_x2 = Sig[1,1] - t(Sig[1,2:k])%*%solve(Sig[2,2])%*%Sig[1,2:k]
  
  ###Calculate the contributions to the log likelihood.
  x2part = 0
  for(i in 1:j){
    temp = dnorm(x2[,i],mean=mu_x2[,i],sd=sqrt(sig2_x2[i,i]),log = TRUE)
    x2part = x2part + temp
  }
  
  #When y1=0:
  ll0 = pnorm(0,mean=mu_y1_x2,sd=sqrt(sig2_y1_x2),log.p = TRUE) + 
    x2part
  
  #When y1>0:
  ll1 = dnorm(y1, mean=mu_y1_x2, sd=sqrt(sig2_y1_x2),log = TRUE) + 
    x2part
  
  #Combine them, based on y1
  ll = ifelse(censored,ll0,ll1)
  
  #Return the negative log-likelihood, since I will be using optim (which minimizes instead of maximizing)
  -sum(ll)
  
}

### calculates the log likelihood for a tobit IV regression
# includes intercept terms

loglik_tobitiv_int <- function(t){
  
  # re-listify
  sig_u = t[grep('sig_u',names(t))]
  sig_v = t[grep('sig_v',names(t))]
  rho = t[grep('rho',names(t))]
  beta1 = t[grep('beta1',names(t))]
  beta2 = t[grep('beta2',names(t))]
  pi1 = t[grep('pi1',names(t))]
  pi2 = t[grep('pi2',names(t))]
  
  j = length(regs$endog)
  l = length(regs$inst)
  
  regStarts = c(grep('subelem1',names(pi1)),length(pi1)+1)
  v = diff(regStarts)
  pi1 = split(pi1, rep(1:length(v),v))
  
  regStarts = c(grep('subelem1',names(pi2)),length(pi2)+1)
  v = diff(regStarts)
  pi2 = split(pi2, rep(1:length(v),v))
  
  # make sig_err
  a = matrix(diag(j)*(unlist(sig_v))^2, nrow = j)
  b = cbind(rho,a)
  Sig_err = rbind(c(sig_u^2,rho),b)
  if(min(eigen(Sig_err)$values)<=0){
    print("negative eigenvalue in Sigma_err")
    return(Inf)
  }
  if((sig_u<=0)){
    print('negative variance value')
    return(Inf)
  }
  
  #Transformation matrix from (eta, u, v) to (y0*, log y1*, x2) (without the mean)
  A = diag(j+1)
  A[1,] = c(1,beta2)
  
  Sig = A%*%Sig_err%*%t(A)
  if(min(eigen(Sig)$values)<=0){
    print("negative eigenvalue in Sigma")
    return(Inf)
  }
  if(sig_u<=0){
    print('negative variance value')
    return(Inf)
  }
  
  # get the means
  censored = y1<=0
  
  # first get the x2 means so they're in the right shape
  mu_x2 = matrix(c(rep(NA,j*length(y1>0))),ncol = j)
  for(i in 1:j){
    formula = regs$endogReg[[i]]
    #mf = model.frame(formula = formula, na.action = NULL)
    mf = model.frame(formula = formula, na.action = NULL)
    m = dim(mf)[2]
    x <- model.matrix(attr(mf, "terms"), data=mf)
    x1temp = as.matrix(x[,1:(m-l)]); ztemp = as.matrix(x[,(m-l+1):m])
    mu_x2[,i] = x1temp%*%as.matrix(pi1[[i]]) + ztemp%*%as.matrix(pi2[[i]])
  }
  
  # then get the means for the y regressions
  formula = regs$formula
  #mf = model.frame(formula = formula, na.action = NULL)
  mf = model.frame(formula = formula)
  m = dim(mf)[2]
  x <- model.matrix(attr(mf,"terms"),data = mf)
  x1 = as.matrix(x[,1:(m-j)]); x2 = as.matrix(x[,(m-j+1):m])
  mu_y1 = x1%*%beta1 + mu_x2%*%beta2
  
  # now get the conditional means
  #Parameters for x2
  k = dim(Sig)[1]
  sig2_x2 = as.matrix(Sig[2,2:k])
  
  #Parameters for y1star given x2
  mu_y1_x2 = mu_y1 + c(t(Sig[1,2:k])%*%solve(Sig[2,2:k])%*%t(x2-mu_x2))
  sig2_y1_x2 = Sig[1,1] - t(Sig[1,2:k])%*%solve(Sig[2,2])%*%Sig[1,2:k]
  
  ###Calculate the contributions to the log likelihood.
  x2part = 0
  for(i in 1:j){
    temp = dnorm(x2[,i],mean=mu_x2[,i],sd=sqrt(sig2_x2[i,i]),log = TRUE)
    x2part = x2part + temp
  }
  
  #When y1=0:
  ll0 = pnorm(0,mean=mu_y1_x2,sd=sqrt(sig2_y1_x2),log.p = TRUE) + 
    x2part
  
  #When y1>0:
  ll1 = dnorm(y1, mean=mu_y1_x2, sd=sqrt(sig2_y1_x2),log = TRUE) + 
    x2part
  
  #Combine them, based on y1
  ll = ifelse(censored,ll0,ll1)
  
  #Return the negative log-likelihood, since I will be using optim (which minimizes instead of maximizing)
  -sum(ll)
  
}

