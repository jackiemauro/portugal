
loglik_lgiv<-function(t){
  ############ preliminaries ##############
  # get everyone named
  pieces = name.pieces(t)
  
  # reconstitute covariance matrix de-cholesky-ify if you're supposed to
  Sig_err = reconstitute.cov(vals=pieces$cov_start
                             ,num=num_endog
                             ,chol=myChol)
  
  # transform covariance matrix with betas and gammas
  Sig = make.covTrans(Sig_err, num_endog, pieces$gamma, pieces$beta)
  if(all(is.na(Sig))){
    print("Bad Sig, skipping")
    return(-Inf)
  }
  
  ################ get means ######################
  # get log of outcome and its censored values
  logy1 = log(outcome)
  censored = outcome<=0
  
  # get unconditional means of x2'sm
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
  if(cant.solve(sig2_x2)){return(-Inf)}
  mu_y0_x2 = mu_y0 + t(Sig[1,(k+1-j):k]%*%solve(sig2_x2)%*%t(endog_mat-mu_x2))
  sig2_y0_x2 = Sig[1,1] - Sig[1,(k+1-j):k]%*%solve(sig2_x2)%*%Sig[(k+1-j):k,1]
  
  #Parameters for log(y1star) given x2
  mu_y1_x2 = mu_y1 + t(Sig[2,(k+1-j):k]%*%solve(sig2_x2)%*%t(endog_mat-mu_x2))
  sig2_y1_x2 = Sig[2,2] - Sig[2,(k+1-j):k]%*%solve(sig2_x2)%*%Sig[(k+1-j):k,2]
  
  #Parameters for y0star given y1star and x2
  if(cant.solve(Sig[2:k,2:k])){return(-Inf)}
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
  if(sig2_y0_y1x2<0){return(-Inf)}
  
  ll1 = pnorm(0,mean=mu_y0_y1x2, sd=sig2_y0_y1x2,log.p=TRUE, lower.tail = FALSE) +
    dnorm(logy1, mean=mu_y1_x2, sd=sqrt(sig2_y1_x2),log = TRUE) + 
    x2part
  
  #Combine them, based on y1
  ll = ifelse(censored,ll0,ll1)
  
  #Return the negative log-likelihood
  print(-sum(ll, na.rm = T))
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

