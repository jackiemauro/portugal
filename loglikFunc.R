# do rho, tau's, sig_v, sig_u come from Sig_err or Sig?
# i think they come from Sig_err, doing that for now, but check 
# you can't put things into this that aren't going to be optimized: name pieces fn has to change

loglik_lgiv<-function(t){
  
  ############ preliminaries ##############
  # get everyone named
  pieces = name.pieces(t)
  
  # de-cholesky-ify if you're supposed to
  if(myChol == T){
    cov_vals = c(1,pieces$cov_start)
    empty = diag(2+num_endog)
    empty[upper.tri(empty,diag = T)]<-cov_vals
    Sig_err = t(empty)%*%empty
    Sig_err = Sig_err/Sig_err[1,1]
  }
  else{
    print("Warning: no cholesky decomposition makes maximization less stable")
    cov_vals = c(1,pieces$cov_start)
    empty = diag(2+num_endog)
    empty[upper.tri(empty,diag = T)]<-cov_vals
    tempty = t(empty)
    tempty[upper.tri(tempty,diag = T)]<-cov_vals
    Sig_err = tempty
  }
  
  # transform covariance matrix with betas and gammas
  Sig = make.covTrans(Sig_err, num_endog, pieces$gamma, pieces$beta)
  
  ################ get means ######################
  # get log of outcome and its censored values
  logy1 = log(outcome)
  censored = outcome<=0
  
  # get unconditional means of x2's
  n = length(pieces$pi)
  mu_x2 = matrix(c(rep(NA,n*length(outcome))),ncol = n)
  for(i in 1:n){
    mu_x2[,i]= ER_mat[[i]]%*%pieces$pi[[i]]
  }
  
  # get unconditional means of y's
  mu_y0 = y_mat%*%pieces$gamma
  mu_y1 = y_mat%*%pieces$beta
  
  # get conditional means
  k = dim(Sig)[1]
  j = num_endog
  sig2_x2 = as.matrix(Sig[(k+1-j):k,(k+1-j):k])
  
  #Parameters for y0star given x2
  mu_y0_x2 = mu_y0 + t(Sig[1,(k+1-j):k]%*%solve(sig2_x2)%*%t(endog_mat-mu_x2))
  sig2_y0_x2 = Sig[1,1] - Sig[1,(k+1-j):k]%*%solve(sig2_x2)%*%Sig[(k+1-j):k,1]
  
  #Parameters for log(y1star) given x2
  mu_y1_x2 = mu_y1 + t(Sig[2,(k+1-j):k]%*%solve(sig2_x2)%*%t(endog_mat-mu_x2))
  sig2_y1_x2 = Sig[2,2] - Sig[2,(k+1-j):k]%*%solve(sig2_x2)%*%Sig[(k+1-j):k,2]
  
  #Parameters for y0star given y1star and x2
  mu_y0_y1x2 = mu_y0 + t(Sig[1,2:k,drop=FALSE]%*%solve(Sig[2:k,2:k])%*%t(cbind(logy1-mu_y1,endog_mat-mu_x2)))
  sig2_y0_y1x2 = Sig[1,1] - Sig[1,2:k,drop=FALSE]%*%solve(Sig[2:k,2:k])%*%Sig[2:k,1,drop=FALSE]
  
  ############# calculate likelihood #############
  x2part = 0
  for(i in 1:j){
    temp = dnorm(endog_mat[,i],mean=mu_x2[,i],sd=sqrt(sig2_x2[i,i]),log = TRUE)
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

name.pieces<-function(t){  
  cov_start = t[1:len_cov]
  beta = t[(len_cov+1):(len_cov+num_betas)]
  gamma = t[(len_cov+num_betas+1):(len_cov+2*num_betas)]
  pis = list()
  k = len_cov+2*num_betas+1
  for(i in 1:num_endog){
    pis[[i]] = t[k:(k+num_pis[[i]]-1)]
    k = k+num_pis[[i]]
  }
  #pis = t[(len_cov+2*num_betas+1):length(t)]
  
  return(list(cov_start = cov_start
              ,beta = beta
              ,gamma = gamma
              ,pi = pis))
}

make.covTrans <- function(a,num_endog,gamma,beta){
  if(is.matrix(a)){
    Sig_err = a
  }
  else{
    Sig_err = matrix(a, ncol = 2+num_endog, byrow = F)
  }
  
  A = diag(2+num_endog)
  A[1,3:(2+num_endog)] <- tail(gamma,num_endog)
  A[2,3:(2+num_endog)] <- tail(beta,num_endog)
  
  Sig = A%*%Sig_err%*%t(A)
  
  if(min(eigen(Sig)$values)<=0){
    print("bad start sigma values")
    stop
  }
  
  return(Sig)
}