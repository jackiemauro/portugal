loglik_cragg <- function(params){
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
  
  ## Get variance
  if(sig_u*sig_v <= 0){return(Inf)}
  Sig_err = matrix(c(1,0,tau0,
                     0,sig_u^2,tau1,
                     tau0,tau1,sig_v^2),
                   ncol = 3, byrow = T)
  if(min(eigen(Sig_err)$values)<=0){return(Inf)}
  
  #Transformation matrix from (eta, u, v) to (y0*, y1*, x2) (without the mean)
  A = rbind(
    c(1,0,gamma2),
    c(0,1,beta2),
    c(0,0,1)
  )
  
  ## Variance matrix
  Sig = A%*%Sig_err%*%t(A)
  
  ## Unconditional means
  mu_x2 = pi1*x1 + pi2*z
  mu_y0 = gamma1*x1 + gamma2*mu_x2
  mu_y1 = beta1*x1 + beta2*mu_x2
  
  ## Variance terms
  sig_x2 = Sig[3,3]
  sig_y0_x2 = Sig[1,1] - Sig[1,3]%*%solve(Sig[3,3])%*%Sig[3,1]
  sig_y1_x2 = Sig[2,2] - Sig[2,3]%*%solve(Sig[3,3])%*%Sig[3,2]
  sig_y0_y1x2 = Sig[1,1] - Sig[1,2:3]%*%solve(Sig[2:3,2:3])%*%Sig[2:3,1]
  
  ## Conditional means
  mu_y0_x2 = mu_y0 + Sig[1,3]%*%solve(Sig[3,3])*(x2 - mu_x2)
  mu_y1_x2 = mu_y1 + Sig[2,3]%*%solve(Sig[3,3])*(x2 - mu_x2)
  mu_y0_y1x2 = mu_y0 + Sig[1,2:3]%*%solve(Sig[2:3,2:3])%*%rbind(y1 - mu_y1,x1-mu_x2)
  
  ## Likelihoods
  ll0 = pnorm(0, mean = mu_y0_x2, sd = sig_y0_x2, log.p = T) +
    dnorm(x2, mean = mu_x2, sd = sig_x2, log = T)
  
  ll1 = pnorm(0, mean = mu_y0_x2, sd = sig_y0_x2, log.p = T, lower.tail = F)+
    dnorm(y1, mean = mu_y1_x2, sd = sig_y1_x2, log = T)+
    dnorm(x2, mean = mu_x2, sd = sig_x2, log = T) - 
    pnorm(0, mean = mu_y1_x2, sd = sig_y1_x2, log.p = T, lower.tail = F)
  
  ll = ifelse(y1<=0,ll0,ll1)
  -sum(ll)
}

loglik_cragg_mat <- function(params){
  
  beta1 = params[grep("beta1",names(params))]
  gamma1 = params[grep("gamma1",names(params))]
  pi1 = params[grep("pi1",names(params))]
  beta2 = params[grep("beta2",names(params))]
  gamma2 = params[grep("gamma2",names(params))] 
  pi2 = params[grep("pi2",names(params))]
  sig_v = params[grep("sig_v",names(params))]
  sig_u = params["sig_u"]
  tau0 = params[grep("tau0",names(params))]
  tau1 = params[grep("tau1",names(params))]
  
  m = length(beta1) #exogenous variables
  j = length(beta2) #endogenous variables
  k = length(pi2) #instrument
  
  if(j>k){
    print("Error: More endogenous variables than instruments")
    stop
  }
  
  if(k>j){print("Warning: Overidentification (more instruments than endogenous variables)")}
  
  
  ## Get variance
  if(min(sig_u,sig_v) <= 0){return(Inf)}
  Sig_x2 = matrix(c(rep(0,j^2)),ncol = j)
  for(ii in 1:j){
    Sig_x2[ii,ii] = sig_v[ii]^2
  }
  
  Sig_y = matrix(c(1,0,0,sig_u^2), ncol = 2, byrow = T)
  Sig_t = matrix(c(tau0,tau1), nrow = 2, byrow = T)
  
  temp = cbind(Sig_y,Sig_t)
  temp2 = cbind(t(Sig_t), Sig_x2)
  Sig_err = rbind(temp,temp2)
  if(min(eigen(Sig_err)$values)<=0){return(Inf)}
  
  #Transformation matrix from (eta, u, v) to (y0*, y1*, x2) (without the mean)
  A = diag(2+j)
  A[1,] = c(1,0,gamma2)
  A[2,] = c(0,1,beta2)
  
  ## Variance matrix
  Sig = A%*%Sig_err%*%t(A)
  
  ## Unconditional means
  mu_x2 = matrix(c(rep(0,j*dim(x1)[2])),ncol = dim(x1)[2])
  for(ii in 1:j){
    mu_x2[ii,] = t(pi1)%*%x1 + pi2[ii]*matrix(z, nrow = k, byrow = T)[ii,]
  }
  mu_y0 = t(gamma1)%*%x1 + t(gamma2)%*%mu_x2
  mu_y1 = t(beta1)%*%x1 + t(beta2)%*%mu_x2
  
  ## Variance terms
  sig2_x2 = Sig[(3:dim(Sig)[1]),(3:dim(Sig)[1])]
  sig2_y0_x2 = Sig[1,1] - Sig[1,(3:dim(Sig)[1])]%*%solve(sig2_x2)%*%Sig[(3:dim(Sig)[1]),1]
  sig2_y1_x2 = Sig[2,2] - Sig[2,(3:dim(Sig)[1])]%*%solve(sig2_x2)%*%Sig[(3:dim(Sig)[1]),2]
  sig2_y0_y1x2 = Sig[1,1] - Sig[1,(2:dim(Sig)[1]),drop=F]%*%solve(Sig[(2:dim(Sig)[1]),(2:dim(Sig)[1])])%*%Sig[(2:dim(Sig)[1]),1,drop=F]
  
  ## Conditional means
  mu_y0_x2 = mu_y0 + Sig[1,(3:dim(Sig)[1])]%*%solve(sig2_x2)%*%matrix(c(x2 - mu_x2),nrow=k, byrow = T)
  mu_y1_x2 = mu_y1 + Sig[2,(3:dim(Sig)[1])]%*%solve(sig2_x2)%*%matrix(c(x2 - mu_x2),nrow=k, byrow = T)
  mu_y0_y1x2 = mu_y0 + Sig[1,(2:dim(Sig)[1]),drop=F]%*%solve(Sig[(2:dim(Sig)[1]),(2:dim(Sig)[1])])%*%rbind(y1 - mu_y1,x2-mu_x2)
  
  ## Likelihoods
  ll0 = pnorm(0, mean = mu_y0_x2, sd = sqrt(sig2_y0_x2), log.p = T) +
    colSums(dnorm(x2, mean = mu_x2, sd = sqrt(sig2_x2), log = T))
  
  ll1 = pnorm(0, mean = mu_y0_x2, sd = sqrt(sig2_y0_x2), log.p = T, lower.tail = F)+
    dnorm(y1, mean = mu_y1_x2, sd = sqrt(sig2_y1_x2), log = T)+
    colSums(dnorm(x2, mean = mu_x2, sd = sqrt(sig2_x2), log = T)) - 
    pnorm(0, mean = mu_y1_x2, sd = sqrt(sig2_y1_x2), log.p = T, lower.tail = F)
  
  ll = ifelse(y1<=0,ll0,ll1)
  ll = ifelse(is.infinite(ll),-1e10,ll)
  -sum(ll)
}