
cragg_errs<-function(cov,df,pi,x1,gamma,beta,n,z){
  
  j=1
  endog = c(rep(NA,n))
  y0 = c(rep(NA,n))
  yStar = c(rep(NA,n))
  errors = matrix(c(rep(0,n*dim(cov)[1])),ncol = dim(cov)[1])
  
  while(j<=n){
    err = mvrnorm(1,rep(0,dim(cov)[1]),cov)
    frame = as.matrix(cbind(1,x1[j,],z[j,]))
    endog[j] = frame%*%pi + err[3]
    frame2 = as.matrix(cbind(1,x1[j,],endog[j]))
    y0[j] = as.numeric(frame2%*%gamma + err[1] > 0)
    yStar[j] = frame2%*%beta + err[2]
    
    if(yStar[j]>0){
      errors[j,] = err
      j = j+1
    }
    
  }
  
  return(list(errors = errors, endog = endog, y0 = y0, yStar = yStar))
}



cragg_errs2<-function(cov,df,pi,x1,gamma,beta,n,z){
  
  newcov = matrix(c(cov[1,1],cov[1,3],cov[1,3],cov[3,3]), ncol = 2, byrow = T)
  err1 = mvrnorm(n,rep(0,2),newcov)
  frame = as.matrix(cbind(1,x1,z))
  endog = frame%*%pi + err1[,2]
  frame2 = as.matrix(cbind(1,x1,endog))
  y0 = as.numeric(frame2%*%gamma + err1[,1] > 0)
  
  # get conditional mean and variance
  mu = c(cov[1,2], cov[2,3])%*%solve(newcov)%*%t(err1)
  sig2 = cov[2,2] - c(cov[1,2], cov[2,3])%*%solve(newcov)%*%c(cov[2,3], cov[1,2])
  
  j=1
  yStar = c(rep(NA,n))
  err2 = c(rep(NA,n))
  
  while(j<=n){
    err = rnorm(1,mu,sqrt(sig2))
    yStar[j] = frame2[j,]%*%beta + err
    if(yStar[j]>0){
      err2[j] = err
      j = j+1
    }
    
  }
  errors = cbind(err1,err2)
  return(list(errors = errors, endog = endog, y0 = y0, yStar = yStar))
}

