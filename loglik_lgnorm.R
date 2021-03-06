### calculates the log likelihood for a lognormal IV regression
# does not include intercept terms

# To do: cholesky decomposition
# better way to find pi's
# need a way to deal with NA's -- leaving them in does terrible things
# replacing with option drop.unused.levels = T does not work
# you get the wrong dimensions

loglik_lgnorm <- function(t){
  
  # re-listify
  sig_u = t[grep('sig_u',names(t))]
  sig_v = t[grep('sig_v',names(t))]
  tau0 = t[grep('tau0',names(t))]
  tau1 = t[grep('tau1',names(t))]
  #rho = t[grep('rho',names(t))]
  beta1 = t[grep('beta1',names(t))]
  beta2 = t[grep('beta2',names(t))]
  gamma1 = t[grep('gamma1',names(t))]
  gamma2 = t[grep('gamma2',names(t))]
  pi1 = t[grep('pi1',names(t))]
  pi2 = t[grep('pi2',names(t))]
  
  j = length(regs$endog)
  l = length(regs$inst)
  logy1 = log(y1)
  
  regStarts = c(grep('subelem1',names(pi1)),length(pi1)+1)
  v = diff(regStarts)
  pi1 = split(pi1, rep(1:length(v),v))
  
  regStarts = c(grep('subelem1',names(pi2)),length(pi2)+1)
  v = diff(regStarts)
  pi2 = split(pi2, rep(1:length(v),v))
  
  #params = t
  # make sig_err
  a = rbind(tau0,tau1)
  #b = matrix(c(1,rho,rho,sig_u),ncol = 2, byrow = T)
  b = matrix(c(1,0,0,sig_u^2),ncol = 2, byrow = T)
  c = matrix(diag(j)*(unlist(sig_v))^2, nrow = j)
  Sig_err = rbind(cbind(b,a),cbind(t(a),c))
  if(min(eigen(Sig_err)$values)<=0){
    print("negative eigenvalue in Sigma_err")
    return(Inf)
  }
  if((sig_u<=0)|(any(sig_v<=0))){
    print('negative variance value')
    return(Inf)
  }
  
  #Transformation matrix from (eta, u, v) to (y0*, log y1*, x2) (without the mean)
  A = diag(j+2)
  A[1,] = c(1,0,gamma2)
  A[2,] = c(0,1,beta2)
  
  Sig = A%*%Sig_err%*%t(A)
  if(min(eigen(Sig)$values)<=0){
    print("negative eigenvalue in Sigma")
    return(Inf)
  }
  if((sig_u<=0)|(any(sig_v<=0))){
    print('negative variance value')
    return(Inf)
  }
  
  # get the means
  mfBig= model.frame(formula = regs$formula)
  drop = which(is.na(rowSums(mfBig)))
  censored = y1[-drop]<=0
  
  # drop all rows with NA in any regression
  mf = list()
  for(i in 1:j){
    formula = regs$endogReg[[i]]
    mf[[i]] = model.frame(formula = formula)
    drop = union(drop,which(is.na(rowSums(mf))))
  }
  
  # first get the x2 means so they're in the right shape
  mu_x2 = matrix(c(rep(NA,j*dim(mfBig[keep,])[1])),ncol = j)
  for(i in 1:j){
    formula = regs$endogReg[[i]]
    mf = model.frame(formula = formula)
    drop = union(drop,which(is.na(rowSums(mf))))
    mf = mf[-drop,]
    m = dim(mf)[2]
    x <- model.matrix(attr(mf, "terms"), data=mf)
    x1temp = as.matrix(x[keep,2:(m-l)])
    ztemp = as.matrix(x[keep,(m-l+1):m])
    mu_x2[,i] = x1temp%*%as.matrix(pi1[[i]]) + ztemp%*%as.matrix(pi2[[i]])
  }
  
  # then get the means for the y regressions
  formula = regs$formula
  mf = model.frame(formula = formula, drop.unused.levels = T)
  mf = mf[keep,]
  m = dim(mf)[2]
  x <- model.matrix(attr(mf,"terms"),data = mf)
  x1 = as.matrix(x[keep,2:(m-j)]); x2 = as.matrix(x[keep,(m-j+1):m])
  mu_y0 = x1%*%gamma1 + mu_x2%*%gamma2
  mu_y1 = x1%*%beta1 + mu_x2%*%beta2
  
  # now get the conditional means
  #Parameters for x2
  k = dim(Sig)[1]
  sig2_x2 = as.matrix(Sig[(k+1-j):k,(k+1-j):k])
  
  #Parameters for y0star given x2
  mu_y0_x2 = mu_y0 + t(Sig[1,(k+1-j):k]%*%solve(sig2_x2)%*%t(x2-mu_x2))
  sig2_y0_x2 = Sig[1,1] - Sig[1,(k+1-j):k]%*%solve(sig2_x2)%*%Sig[(k+1-j):k,1]
  
  #Parameters for log(y1star) given x2
  mu_y1_x2 = mu_y1 + t(Sig[2,(k+1-j):k]%*%solve(sig2_x2)%*%t(x2-mu_x2))
  sig2_y1_x2 = Sig[2,2] - Sig[2,(k+1-j):k]%*%solve(sig2_x2)%*%Sig[(k+1-j):k,2]
  
  #Parameters for y0star given y1star and x2
  mu_y0_y1x2 = mu_y0 + t(Sig[1,2:k,drop=FALSE]%*%solve(Sig[2:k,2:k])%*%t(cbind(logy1[keep]-mu_y1,x2-mu_x2)))
  sig2_y0_y1x2 = Sig[1,1] - Sig[1,2:k,drop=FALSE]%*%solve(Sig[2:k,2:k])%*%Sig[2:k,1,drop=FALSE]
  
  ###Calculate the contributions to the log likelihood.
  x2part = 0
  for(i in 1:j){
    temp = dnorm(x2[,i],mean=mu_x2[,i],sd=sqrt(sig2_x2[i,i]),log = TRUE)
    x2part = x2part + temp
  }
  #When y1=0:
  ll0 = pnorm(0,mean=mu_y0_x2,sd=sqrt(sig2_y0_x2),log.p = TRUE) + x2part
  
  #When y1>0:
  ll1 = pnorm(0,mean=mu_y0_y1x2, sd=sqrt(sig2_y0_y1x2),log.p=TRUE, lower.tail = FALSE) +
    dnorm(logy1[keep], mean=mu_y1_x2, sd=sqrt(sig2_y1_x2),log = TRUE) + 
    x2part
  
  #Combine them, based on y1
  ll = ifelse(censored,ll0,ll1)
  
  #Return the negative log-likelihood
  -sum(ll, na.rm = T)
  
}

### calculates the log likelihood for a lognormal IV regression
# does not include intercept terms

# To do: cholesky decomposition
# better way to find pi's

loglik_lgnorm_int <- function(t){
  
  # re-listify
  sig_u = t[grep('sig_u',names(t))]
  sig_v = t[grep('sig_v',names(t))]
  tau0 = t[grep('tau0',names(t))]
  tau1 = t[grep('tau1',names(t))]
  #rho = t[grep('rho',names(t))]
  beta1 = t[grep('beta1',names(t))]
  beta2 = t[grep('beta2',names(t))]
  gamma1 = t[grep('gamma1',names(t))]
  gamma2 = t[grep('gamma2',names(t))]
  pi1 = t[grep('pi1',names(t))]
  pi2 = t[grep('pi2',names(t))]
  
  j = length(regs$endog)
  l = length(regs$inst)
  logy1 = log(y1)
  
  regStarts = c(grep('subelem1',names(pi1)),length(pi1)+1)
  v = diff(regStarts)
  pi1 = split(pi1, rep(1:length(v),v))
  
  regStarts = c(grep('subelem1',names(pi2)),length(pi2)+1)
  v = diff(regStarts)
  pi2 = split(pi2, rep(1:length(v),v))
  
  #params = t
  # make sig_err
  a = rbind(tau0,tau1)
  #b = matrix(c(1,rho,rho,sig_u),ncol = 2, byrow = T)
  b = matrix(c(1,0,0,sig_u^2),ncol = 2, byrow = T)
  c = matrix(diag(j)*(unlist(sig_v))^2, nrow = j)
  Sig_err = rbind(cbind(b,a),cbind(t(a),c))
  if(min(eigen(Sig_err)$values)<=0){
    print("negative eigenvalue in Sigma_err")
    return(Inf)
  }
  if((sig_u<=0)|(any(sig_v<=0))){
    print('negative variance value')
    return(Inf)
  }
  
  #Transformation matrix from (eta, u, v) to (y0*, log y1*, x2) (without the mean)
  A = diag(j+2)
  A[1,] = c(1,0,gamma2)
  A[2,] = c(0,1,beta2)
  
  Sig = A%*%Sig_err%*%t(A)
  if(min(eigen(Sig)$values)<=0){
    print("negative eigenvalue in Sigma")
    return(Inf)
  }
  if((sig_u<=0)|(any(sig_v<=0))){
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
    mf = model.frame(formula = formula, drop.unused.levels =T)
    m = dim(mf)[2]
    x <- model.matrix(attr(mf, "terms"), data=mf)
    x1temp = as.matrix(x[,1:(m-l)]); ztemp = as.matrix(x[,(m-l+1):m])
    mu_x2[,i] = x1temp%*%as.matrix(pi1[[i]]) + ztemp%*%as.matrix(pi2[[i]])
  }
  
  # then get the means for the y regressions
  formula = regs$formula
  #mf = model.frame(formula = formula, na.action = NULL)
  mf = model.frame(formula = formula, drop.unused.levels = T)
  m = dim(mf)[2]
  x <- model.matrix(attr(mf,"terms"),data = mf)
  x1 = as.matrix(x[,1:(m-j)]); x2 = as.matrix(x[,(m-j+1):m])
  mu_y0 = x1%*%gamma1 + mu_x2%*%gamma2
  mu_y1 = x1%*%beta1 + mu_x2%*%beta2
  
  # now get the conditional means
  #Parameters for x2
  k = dim(Sig)[1]
  sig2_x2 = as.matrix(Sig[(k+1-j):k,(k+1-j):k])
  
  #Parameters for y0star given x2
  mu_y0_x2 = mu_y0 + t(Sig[1,(k+1-j):k]%*%solve(sig2_x2)%*%t(x2-mu_x2))
  sig2_y0_x2 = Sig[1,1] - Sig[1,(k+1-j):k]%*%solve(sig2_x2)%*%Sig[(k+1-j):k,1]
  
  #Parameters for log(y1star) given x2
  mu_y1_x2 = mu_y1 + t(Sig[2,(k+1-j):k]%*%solve(sig2_x2)%*%t(x2-mu_x2))
  sig2_y1_x2 = Sig[2,2] - Sig[2,(k+1-j):k]%*%solve(sig2_x2)%*%Sig[(k+1-j):k,2]
  
  #Parameters for y0star given y1star and x2
  mu_y0_y1x2 = mu_y0 + t(Sig[1,2:k,drop=FALSE]%*%solve(Sig[2:k,2:k])%*%t(cbind(logy1-mu_y1,x2-mu_x2)))
  sig2_y0_y1x2 = Sig[1,1] - Sig[1,2:k,drop=FALSE]%*%solve(Sig[2:k,2:k])%*%Sig[2:k,1,drop=FALSE]
  
  ###Calculate the contributions to the log likelihood.
  x2part = 0
  for(i in 1:j){
    temp = dnorm(x2[,i],mean=mu_x2[,i],sd=sqrt(sig2_x2[i,i]),log = TRUE)
    x2part = x2part + temp
  }
  #When y1=0:
  ll0 = pnorm(0,mean=mu_y0_x2,sd=sqrt(sig2_y0_x2),log.p = TRUE) + x2part
  
  #When y1>0:
  ll1 = pnorm(0,mean=mu_y0_y1x2, sd=sqrt(sig2_y0_y1x2),log.p=TRUE, lower.tail = FALSE) +
    dnorm(logy1, mean=mu_y1_x2, sd=sqrt(sig2_y1_x2),log = TRUE) + 
    x2part
  
  #Combine them, based on y1
  ll = ifelse(censored,ll0,ll1)
  
  #Return the negative log-likelihood
  -sum(ll, na.rm = T)
  
}
