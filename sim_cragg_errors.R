
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
    y0[j] = as.numeric(frame%*%gamma + err[1] > 0)
    yStar[j] = frame2%*%beta + err[2]
    
    if(yStar[j]>0){
      errors[j,] = err
      j = j+1
    }
    
  }
  
  return(list(errors = errors, endog = endog, y0 = y0, yStar = yStar))
}
