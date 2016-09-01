### This function returns the negative log-likelihood (so that it can be minimized)
### The math follows the derivation in the pdf I sent over
### There is not much error checking on the parameters, so be careful.
### The function does check for a PSD (eta,u,v) covariance matrix, and returns Inf if it fails.  This works well with optim.

## still to do: program in relaxation of independence
loglik_lgnormO = function(params){
  sig_v = params[1]
  sig_u = params[2]
  tau0 = params[3]
  tau1 = params[4]
  beta1 = params[5]
  beta2 = params[6]
  gamma1 = params[7]
  gamma2 = params[8]
  pi1 = params[9]
  pi2 = params[10]
  
  censored = y1<=0
  #We will just work with the log of y1, so that we don't need to care about lognormal distributions.
  logy1 = log(y1)
  
  #Check for invalid parameters
  Sig_err= matrix( c(1, 0,    tau0,
                     0,     sig_u^2,    tau1,
                     tau0,  tau1, sig_v^2),
                   ncol = 3, byrow = T)
  if(min(eigen(Sig_err)$values)<=0){return(Inf)}
  if((sig_u<=0)|(sig_v<=0)){return(Inf)}
  
  #Transformation matrix from (eta, u, v) to (y0*, log y1*, x2) (without the mean)
  A = rbind(
    c(1,0,gamma2),
    c(0,1,beta2),
    c(0,0,1)
  )
  #Covariance of (y0*, log y1*, x2)
  Sig = A%*%Sig_err%*%t(A)
  
  #### WARNING: I WILL BE SLOPPY FROM HERE ON BY DROPPING THE LOG IN MY log y1* NOTATION (and the star) ####
  #### It's just too much to carry around
  
  #Means for (y0^*, log y1*, x2).
  mu_y0 = (gamma1+pi1*gamma2)*x1 + gamma2*pi2*z
  mu_y1 = (beta1+pi1*beta2)*x1+pi2*beta2*z
  mu_x2 = x1*pi1 + z*pi2
  
  #### NOTATION: For mu and sig2 variables, I will separate the variables from the conditioning variables by underscores.
  #### That is, (mu/sig2)_(Variable)_(Conditioning variables)
  
  #Parameters for x2
  sig2_x2 = Sig[3,3]
  
  #Parameters for y0star given x2
  mu_y0_x2 = mu_y0 + Sig[1,3]/Sig[3,3]*(x2-mu_x2)
  sig2_y0_x2 = Sig[1,1] - Sig[1,3]^2/Sig[3,3]
  
  #Parameters for log(y1star) given x2
  mu_y1_x2 = mu_y1 + Sig[2,3]/Sig[3,3]*(x2-mu_x2)
  sig2_y1_x2 = Sig[2,2] - Sig[2,3]^2/Sig[3,3]
  
  #Parameters for y0star given y1star and x2
  mu_y0_y1x2 = mu_y0 + Sig[1,2:3,drop=FALSE]%*%solve(Sig[2:3,2:3])%*%rbind(logy1-mu_y1,x2-mu_x2)
  sig2_y0_y1x2 = Sig[1,1] - Sig[1,2:3,drop=FALSE]%*%solve(Sig[2:3,2:3])%*%Sig[2:3,1,drop=FALSE]
  
  ###Calculate the contributions to the log likelihood.
  
  #When y1=0:
  ll0 = pnorm(0,mean=mu_y0_x2,sd=sqrt(sig2_y0_x2),log.p = TRUE) + 
    dnorm(x2,mean=mu_x2,sd=sqrt(sig2_x2),log = TRUE)
  
  #When y1>0:
  ll1 = pnorm(0,mean=mu_y0_y1x2, sd=sqrt(sig2_y0_y1x2),log.p=TRUE, lower.tail = FALSE) +
    dnorm(logy1, mean=mu_y1_x2, sd=sqrt(sig2_y1_x2),log = TRUE) + 
    dnorm(x2,mean=mu_x2,sd=sqrt(sig2_x2),log = TRUE)
  
  #Combine them, based on y1
  ll = ifelse(censored,ll0,ll1)
  
  #Return the negative log-likelihood, since I will be using optim (which minimizes instead of maximizing)
  -sum(ll)
}

#eventually you should have this one be the primary function
#keeping both until this is verified to work well
loglik_lgnorm_mat = function(params){
  m = dim(x1)[1]
  sig_v = params[1]
  sig_u = params[2]
  rho = params[3]
  tau0 = params[4]
  tau1 = params[5]
  beta1 = params[6:(5+m)]
  beta2 = params[(6+m)]
  gamma1 = params[(7+m):(6+2*m)]
  gamma2 = params[(7+2*m)]
  pi1 = params[(8+2*m):(7+3*m)]
  pi2 = params[8+3*m]
  
  censored = y1<=0
  #We will just work with the log of y1, so that we don't need to care about lognormal distributions.
  logy1 = log(y1)
  
  #Check for invalid parameters
  Sig_err= matrix( c(1, rho,    tau0,
                     rho,     sig_u^2,    tau1,
                     tau0,  tau1, sig_v^2),
                   ncol = 3, byrow = T)
  if(min(eigen(Sig_err)$values)<=0){return(Inf)}
  if((sig_u<=0)|(sig_v<=0)){return(Inf)}
  
  #Transformation matrix from (eta, u, v) to (y0*, log y1*, x2) (without the mean)
  A = rbind(
    c(1,0,gamma2),
    c(0,1,beta2),
    c(0,0,1)
  )
  #Covariance of (y0*, log y1*, x2)
  Sig = A%*%Sig_err%*%t(A)
  
  #### WARNING: I WILL BE SLOPPY FROM HERE ON BY DROPPING THE LOG IN MY log y1* NOTATION (and the star) ####
  #### It's just too much to carry around
  
  #Means for (y0^*, log y1*, x2).
  mu_y0 = t(gamma1+pi1*gamma2)%*%x1 + gamma2*pi2*z
  mu_y1 = t(beta1+pi1*beta2)%*%x1+pi2*beta2*z
  mu_x2 = t(pi1)%*%x1 + z*pi2
  
  #### NOTATION: For mu and sig2 variables, I will separate the variables from the conditioning variables by underscores.
  #### That is, (mu/sig2)_(Variable)_(Conditioning variables)
  
  #Parameters for x2
  sig2_x2 = Sig[3,3]
  
  #Parameters for y0star given x2
  mu_y0_x2 = mu_y0 + Sig[1,3]/Sig[3,3]*(x2-mu_x2)
  sig2_y0_x2 = Sig[1,1] - Sig[1,3]^2/Sig[3,3]
  
  #Parameters for log(y1star) given x2
  mu_y1_x2 = mu_y1 + Sig[2,3]/Sig[3,3]*(x2-mu_x2)
  sig2_y1_x2 = Sig[2,2] - Sig[2,3]^2/Sig[3,3]
  
  #Parameters for y0star given y1star and x2
  mu_y0_y1x2 = mu_y0 + Sig[1,2:3,drop=FALSE]%*%solve(Sig[2:3,2:3])%*%rbind(logy1-mu_y1,x2-mu_x2)
  sig2_y0_y1x2 = Sig[1,1] - Sig[1,2:3,drop=FALSE]%*%solve(Sig[2:3,2:3])%*%Sig[2:3,1,drop=FALSE]
  
  ###Calculate the contributions to the log likelihood.
  
  #When y1=0:
  ll0 = pnorm(0,mean=mu_y0_x2,sd=sqrt(sig2_y0_x2),log.p = TRUE) + 
    dnorm(x2,mean=mu_x2,sd=sqrt(sig2_x2),log = TRUE)
  
  #When y1>0:
  ll1 = pnorm(0,mean=mu_y0_y1x2, sd=sqrt(sig2_y0_y1x2),log.p=TRUE, lower.tail = FALSE) +
    dnorm(logy1, mean=mu_y1_x2, sd=sqrt(sig2_y1_x2),log = TRUE) + 
    dnorm(x2,mean=mu_x2,sd=sqrt(sig2_x2),log = TRUE)
  
  #Combine them, based on y1
  ll = ifelse(censored,ll0,ll1)
  
  #Return the negative log-likelihood, since I will be using optim (which minimizes instead of maximizing)
  -sum(ll)
}

#keeping both until this is verified to work well
loglik_lgnorm_norho = function(params){
  m = dim(x1)[1]
  sig_v = params[1]
  sig_u = params[2]
  tau0 = params[3]
  tau1 = params[4]
  beta1 = params[5:(4+m)]
  beta2 = params[(5+m)]
  gamma1 = params[(6+m):(5+2*m)]
  gamma2 = params[(6+2*m)]
  pi1 = params[(7+2*m):(6+3*m)]
  pi2 = params[7+3*m]
  
  censored = y1<=0
  #We will just work with the log of y1, so that we don't need to care about lognormal distributions.
  logy1 = log(y1)
  
  #Check for invalid parameters
  Sig_err= matrix( c(1, 0,    tau0,
                     0,     sig_u^2,    tau1,
                     tau0,  tau1, sig_v^2),
                   ncol = 3, byrow = T)
  if(min(eigen(Sig_err)$values)<=0){return(Inf)}
  if((sig_u<=0)|(sig_v<=0)){return(Inf)}
  
  #Transformation matrix from (eta, u, v) to (y0*, log y1*, x2) (without the mean)
  A = rbind(
    c(1,0,gamma2),
    c(0,1,beta2),
    c(0,0,1)
  )
  #Covariance of (y0*, log y1*, x2)
  Sig = A%*%Sig_err%*%t(A)
  
  #Means for (y0*, log y1*, x2).
  mu_y0 = t(gamma1+pi1*gamma2)%*%x1 + gamma2*pi2*z
  mu_y1 = t(beta1+pi1*beta2)%*%x1+pi2*beta2*z
  mu_x2 = t(pi1)%*%x1 + z*pi2
  
  #### NOTATION: For mu and sig2 variables, separate the variables from the conditioning variables by underscores.
  #### That is, (mu/sig2)_(Variable)_(Conditioning variables)
  
  #Parameters for x2
  sig2_x2 = Sig[3,3]
  
  #Parameters for y0star given x2
  mu_y0_x2 = mu_y0 + Sig[1,3]/Sig[3,3]*(x2-mu_x2)
  sig2_y0_x2 = Sig[1,1] - Sig[1,3]^2/Sig[3,3]
  
  #Parameters for log(y1star) given x2
  mu_y1_x2 = mu_y1 + Sig[2,3]/Sig[3,3]*(x2-mu_x2)
  sig2_y1_x2 = Sig[2,2] - Sig[2,3]^2/Sig[3,3]
  
  #Parameters for y0star given y1star and x2
  mu_y0_y1x2 = mu_y0 + Sig[1,2:3,drop=FALSE]%*%solve(Sig[2:3,2:3])%*%rbind(logy1-mu_y1,x2-mu_x2)
  sig2_y0_y1x2 = Sig[1,1] - Sig[1,2:3,drop=FALSE]%*%solve(Sig[2:3,2:3])%*%Sig[2:3,1,drop=FALSE]
  
  ###Calculate the contributions to the log likelihood.
  
  #When y1=0:
  ll0 = pnorm(0,mean=mu_y0_x2,sd=sqrt(sig2_y0_x2),log.p = TRUE) + 
    dnorm(x2,mean=mu_x2,sd=sqrt(sig2_x2),log = TRUE)
  
  #When y1>0:
  ll1 = pnorm(0,mean=mu_y0_y1x2, sd=sqrt(sig2_y0_y1x2),log.p=TRUE, lower.tail = FALSE) +
    dnorm(logy1, mean=mu_y1_x2, sd=sqrt(sig2_y1_x2),log = TRUE) + 
    dnorm(x2,mean=mu_x2,sd=sqrt(sig2_x2),log = TRUE)
  
  #Combine them, based on y1
  ll = ifelse(censored,ll0,ll1)
  
  #Return the negative log-likelihood
  -sum(ll)
}

#Adapting norho version to include multiple instruments
loglik_lgnorm_norho = function(params){
  sig_v = un[1:(breaks[1]-1)]
  sig_u = un[(breaks[1]+1):(breaks[2]-1)]
  tau0 = un[(breaks[2]+1):(breaks[3]-1)]
  tau1 = un[(breaks[3]+1):(breaks[4]-1)]
  beta1 = un[(breaks[4]+1):(breaks[5]-1)]
  beta2 = un[(breaks[5]+1):(breaks[6]-1)]
  gamma1 = un[(breaks[6]+1):(breaks[7]-1)]
  gamma2 = un[(breaks[7]+1):(breaks[8]-1)]
  pi1 = un[(breaks[8]+1):(breaks[9]-1)]
  pi2 = un[(breaks[9]+1):(breaks[10]-1)]
  
  j = length(pi2) # change this to the number of x2's
  m = dim(x1)[1]
  
  censored = y1<=0
  logy1 = log(y1)
  
  #Get unconditional covariance matrix
  #assume x2's (endogenous variables) are uncorrelated
  Sigx2 = diag(j)*sig_v
  int = matrix(c(1,0,0,sig_u^2),ncol=2,byrow=T)
  tauMat = rbind(tau0,tau1)
  Sig_err = rbind(cbind(int,tauMat),cbind(t(tauMat),Sigx2))
  if(min(eigen(Sig_err)$values)<=0){return(Inf)}
  if((sig_u<=0)|(sig_v<=0)){return(Inf)}
  
  #Transformation matrix from (eta, u, v) to (y0*, log y1*, x2) (without the mean)
  A = diag(j+2)
  A[1,] = c(1,0,gamma2)
  A[2,] = c(0,1,beta2)
  
  #Covariance of (y0*, log y1*, x2)
  Sig = A%*%Sig_err%*%t(A)
  
  #Means for (y0*, log y1*, x2)
  mu_x2 = matrix(c(rep(NA,j)), ncol = 1)
  for(i in 1:j){mu_x2[i,] = pi1%*%x1 + pi2[i]*z[i,]}
  mu_y0 = t(gamma1)%*%x1 + t(gamma2)%*%(t(pi1)%*%x1 + t(pi2)%*%z)
  mu_y0 = t(gamma1+pi1*gamma2)%*%x1 + gamma2*pi2*z
  mu_y1 = t(beta1+pi1*beta2)%*%x1+pi2*beta2*z
  
  
  #### NOTATION: For mu and sig2 variables, separate the variables from the conditioning variables by underscores.
  #### That is, (mu/sig2)_(Variable)_(Conditioning variables)
  
  #Parameters for x2
  sig2_x2 = Sig[3,3]
  
  #Parameters for y0star given x2
  mu_y0_x2 = mu_y0 + Sig[1,3]/Sig[3,3]*(x2-mu_x2)
  sig2_y0_x2 = Sig[1,1] - Sig[1,3]^2/Sig[3,3]
  
  #Parameters for log(y1star) given x2
  mu_y1_x2 = mu_y1 + Sig[2,3]/Sig[3,3]*(x2-mu_x2)
  sig2_y1_x2 = Sig[2,2] - Sig[2,3]^2/Sig[3,3]
  
  #Parameters for y0star given y1star and x2
  mu_y0_y1x2 = mu_y0 + Sig[1,2:3,drop=FALSE]%*%solve(Sig[2:3,2:3])%*%rbind(logy1-mu_y1,x2-mu_x2)
  sig2_y0_y1x2 = Sig[1,1] - Sig[1,2:3,drop=FALSE]%*%solve(Sig[2:3,2:3])%*%Sig[2:3,1,drop=FALSE]
  
  ###Calculate the contributions to the log likelihood.
  
  #When y1=0:
  ll0 = pnorm(0,mean=mu_y0_x2,sd=sqrt(sig2_y0_x2),log.p = TRUE) + 
    dnorm(x2,mean=mu_x2,sd=sqrt(sig2_x2),log = TRUE)
  
  #When y1>0:
  ll1 = pnorm(0,mean=mu_y0_y1x2, sd=sqrt(sig2_y0_y1x2),log.p=TRUE, lower.tail = FALSE) +
    dnorm(logy1, mean=mu_y1_x2, sd=sqrt(sig2_y1_x2),log = TRUE) + 
    dnorm(x2,mean=mu_x2,sd=sqrt(sig2_x2),log = TRUE)
  
  #Combine them, based on y1
  ll = ifelse(censored,ll0,ll1)
  
  #Return the negative log-likelihood
  -sum(ll)
}
