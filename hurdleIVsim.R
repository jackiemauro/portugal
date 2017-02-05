

# still needs to be able to incorporate more than one endogenous variable
# still needs to be able to specify which exogenous variables to include when


hurdle.IV.sim <- function(formula = F,
                          pi = c(1,-1,3), #intercept, exog, inst
                          gamma = c(-.2,.8,.07), #intercept, exog, endog
                          beta = c(.05,.06,.02), #intercept, exog, endog
                          endog_reg = list(),
                          exog_mean = 1,
                          exog_sd = 1,
                          z_mean = 3,
                          z_sd = 1,
                          endog_sd = 5,
                          y_sd = 2,
                          rho = .2,
                          tau0 = .3,
                          tau1 = .1,
                          n = 1000,
                          silent = F,
                          type = "lognormal"){
  
  #checks
  n_z = length(z_mean)
  n_x2 = length(endog_sd)
  
  if(n_x2 > n_z){
    print("Error: more endogenous variables than instruments")
    stop
  }
  if(length(gamma)!= length(beta) ){
    print("Error: coefficient vector lengths differ")
    stop
  }
  if(length(tau0)!=n_x2|length(tau1)!=n_x2){
    print("Error: length of tau vectors must equal number of endogenous variables")
    stop
  }
  
  # make instruments and exogenous variables as random normals
  z = make.df(mean = z_mean, sd = z_sd, n = n, pref = "inst")
  x1 = make.df(mean = exog_mean, sd = exog_sd, n = n, pref = "exog")
  df = data.frame(x1,z)
  
  # make covariance matrix
  require(MASS)
  cov = make.covTrans(list(rho = rho, y_sd = y_sd,endog_sd = endog_sd,
                           tau0 = tau0,tau1 = tau1)
                      , num_endog = length(endog_sd)
                      , gamma = gamma, beta = beta
                      , option = "parameters"
                      , noname = T)
  errors = mvrnorm(n,rep(0,dim(cov)[1]),cov)
  
  # make endogenous variables
  if(length(endog_reg)!=0){
    # need to allow you to exclude some exogenous variables
    print("This functionality not developed yet, set endog_reg = list()")
    stop
  }
  else{
    frame = as.matrix(cbind(1,df))
    endog = frame%*%pi + errors[,3]
  }
  
  # make y
  if(formula!=F){
    # need to allow you to exclude some exogenous variables
    print("This functionality not developed yet, set formula = F")
    stop
  }
  else{
    frame = as.matrix(cbind(1,x1,endog))
    y0 = as.numeric(frame%*%gamma + errors[,1] > 0)
    yStar = frame%*%beta + errors[,2]
    if(type == "lognormal"){
      y = exp(yStar)*y0
    }
    else if(type == "cragg"){
      # should this have as sd y_sd or sd(errors[,2])?
      require(truncnorm)
      y1t = rtruncnorm(n,a=0,mean = frame%*%beta, sd = y_sd)
      y = y1t*y0
    }
    
  }
  out = data.frame(y,y0,endog,df)
  
  if(silent == F){
    par(mfrow = c(1,2))
    hist(log(out$y), main = "log(y)",xlab = "")
    hist(out$y, main = "y", xlab = "")
    par(mfrow = c(1,1))
    print(paste(mean(out$y0)*100," percent uncensored"))
  }
  return(out)
}

make.df <- function(mean,sd,n, pref = F){
  if(length(mean)!=length(sd)){
    print("Error: length of mean and sd vectors differ")
    return(NA)
  }
  out= matrix(c(rep(NA,length(mean)*n)),ncol = length(mean))
  for(i in 1:length(mean)){
    out[,i] = rnorm(n,mean[i],sd[i])
  }
  
  df = as.data.frame(out)
  
  if(pref != F){
    names(df) <- make.names(len = length(mean),pref = pref)
  }
  
  return(df)
}

make.names <- function(len, pref){
  name = c(rep(NA,len))
  for(i in 1:len){
    name[i] = paste(pref,i,sep = "")
  }
  return(name)
}

make.covTrans1 <- function(rho,tau0,tau1,y_sd,endog_sd,gamma,beta){
  mat1 = matrix(c(1,rho,rho,y_sd^2),ncol = 2, byrow = T)
  tau_mat = as.matrix(cbind(tau0,tau1))
  endog_mat = diag(length(endog_sd))*endog_sd^2
  Sig_err = rbind(cbind(mat1,t(tau_mat)),cbind(tau_mat,endog_mat))
  
  m = length(endog_sd)
  A = diag(2+m)
  A[1,3:(2+m)] <- tail(gamma,m)
  A[2,3:(2+m)] <- tail(beta,m)
  
  Sig = A%*%Sig_err%*%t(A)
  
  if(min(eigen(Sig)$values)<=0){
    print("bad start sigma values")
    stop
  }
  
  return(Sig)
}
