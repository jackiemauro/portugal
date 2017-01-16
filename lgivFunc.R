# get start stuff into a vector for the loglikelihood function

lgiv<-function(formula
               , endog_reg
               , data
               , start
               , cholesky = T
               , run_until_PSD = T
               , maxit = 5000
               , trace = 1
               , method = "BFGS"){
  
  #### get covariance matrix ####
  mat1 = matrix(c(1,start$rho,start$rho,start$y_sd^2),ncol = 2, byrow = T)
  tau_mat = as.matrix(cbind(start$tau0,start$tau1))
  endog_mat = diag(length(start$endog_sd))*start$endog_sd^2
  start_cov = rbind(cbind(mat1,t(tau_mat)),cbind(tau_mat,endog_mat))
  
  # transform covariance matrix with beta's and gamma's
  trans_cov = make.covTrans(a = start_cov
                            , num_endog = length(start$endog_sd)
                            , beta = start$beta
                            , gamma = start$gamma)
  
  if(cholesky == T){
    trans_cov = chol(trans_cov)
  }
  else{
    print("Warning: without Cholesky decomposition, maximization is less stable")
    trans_cov = c(trans_cov)
  }
  
  start_vec = c(length(start$beta)
                , length(start$gamma)
                , length(start$endog_sd)
                , )
  out = optim(start_vec
              , loglik_lgiv
              , method= method
              , hessian  = T
              , control = list(maxit = maxit,trace = trace)
              )
  

}