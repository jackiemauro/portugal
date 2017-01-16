# do rho, tau's, sig_v, sig_u come from Sig_err or Sig?
# i think they come from Sig_err, doing that for now, but check 

loglik_lgiv<-function(t){
  
  # get everyone named
  pieces = name.pieces(t)
  
  # cholesky-ify if you're supposed to
  if(pieces$chol == T){
    Sig_err = t(pieces$cov_start)%*%pieces$cov_start
    Sig_err = Sig_err/Sig_err[1,1]
  }
  else{
    print("Warning: no cholesky decomposition makes maximization less stable")
    Sig_err = pieces$cov_start
  }
  
  # transform with betas and gammas
  Sig = make.covTrans(Sig_err, pieces$num_endogs, pieces$gamma, pieces$beta)
  
}

name.pieces<-function(t){
  num_betas = t[1]
  num_gammas = t[1]
  num_endogs = t[2]
  len_cov = (2+num_endogs)^2
  num_pis = NULL
  for(i in 1:num_endogs){
    num_pis[i] = t[2+i]
  }
  k = 2+num_endogs+len_cov
  cov_start = matrix(t[(3+num_endogs):k]
                     , ncol = 2+num_endogs
                     , byrow = F)
  beta = t[(k+1):(k+num_betas)]
  gamma = t[(k+num_betas+1):(k+2*num_betas)]
  pis = t[(k+2*num_betas+1):(length(t)-1)]
  chol = ifelse(tail(t,1)==1,T,F)
  
  return(list(cov_start = cov_start
              ,beta = beta
              ,gamma = gamma
              ,pi = pis
              ,chol = chol
              ,num_endogs = num_endogs))
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