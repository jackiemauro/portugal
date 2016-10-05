#################################################################
# Combining data simulation code into one
# putting too much variance into the x2's
#################################################################


hurdleIV.gen_hurdleSim <- function(formula,
                          family,
                          params,
                          het = FALSE,
                          clump = FALSE,
                          n=10000,
                          rho=F,
                          cond = T){
  
  require(MASS)
  
  mux1 = params$mux1; sigx1 = params$sigx1
  muz = params$muz; sigz = params$sigz
  beta1 = params$beta1; gamma1 = params$gamma1; pi1 = params$pi1
  beta2 = params$beta2; gamma2 = params$gamma2; pi2 = params$pi2
  sig_v = params$sig_v; sig_u = params$sig_u
  tau0 = params$tau0; tau1 = params$tau1
  m = length(formula$exog) #number of exogenous variables
  j = length(formula$inst) #number of instruments
  k = length(formula$endog) #number of endogenous variables
  
  if(j<k){
    print("Error: Too few instruments")
    return(NA)
  }
  
  if(j>k){
    print("Warning: More instruments than endogenous variables")
  }
  
  if(het == TRUE){
    # FIGURE OUT HOW TO PROGRAM IN HETEROSCEDASTICITY
  }
  
  Sig_x2 = matrix(diag(k)*(unlist(sig_v))^2, nrow = k)
  
  if(rho==F){rho=0}
  Sig_y = matrix(c(1,rho,rho,sig_u^2), ncol = 2, byrow = T)
  Sig_t = matrix(c(tau0,tau1), nrow = 2, byrow = T)
  
  temp = cbind(Sig_y,Sig_t)
  temp2 = cbind(t(Sig_t), Sig_x2)
  Sig_err = rbind(temp,temp2)
  
  if(min(eigen(Sig_err)$values)<=0){
    print("bad start sigma values")
    stop
  }
  if(min(sig_v,sig_u)<=0){
    print("bad start sigma values")
    stop
  }
  
  if(cond == F){
    euv = mvrnorm(n=n, mu = c(rep(0,dim(Sig_err)[1])), Sigma = Sig_err)
    eta = euv[,1]
    u = euv[,2]
    v = matrix(c(euv[,3:dim(euv)[2]]),ncol = k)
  }
  else{
    #Transformation matrix from (eta, u, v) to (y0*, log y1*, x2) (without the mean)
    A = diag(k+2)
    A[1,] = c(1,0,beta2)
    A[2,] = c(0,1,gamma2)
    
    Sig = A%*%Sig_err%*%t(A)
    
    #Generate a matrix of the eta,u,v errors and split them out
    euv = mvrnorm(n=n, mu = c(rep(0,dim(Sig)[1])), Sigma = Sig)
    eta = euv[,1]
    u = euv[,2]
    v = matrix(c(euv[,3:dim(euv)[2]]),ncol = k)
  }

  #Construct the variables
  x1 = matrix(c(rep(NA,m*n)),ncol=m)
  nams = c(rep(NA,m))
  for(ii in 1:m){
    nams[ii] <- paste("x1", ii, sep = "")
    x1[,ii] = c(rnorm(n,mux1[(ii)],sigx1[(ii)]))
    assign(nams[ii], x1[,ii])
  }
  colnames(x1) = nams
  
  z = matrix(c(rep(NA,j*n)),ncol=j)
  nams = c(rep(NA,j))
  for(ii in 1:j){
    nams[ii] = paste("z", ii, sep = "")
    z[,ii] = c(rnorm(n,muz[ii],sigz[ii]))
    assign(nams[ii], z[,ii])
  }
  colnames(z) = nams
  
  
  x2 = matrix(c(rep(NA,k*n)),ncol=k)
  nams = c(rep(NA,k))
  for(ii in 1:k){
    nams[ii] = paste("x2", ii, sep = "")
    assign(nams[ii], c(rep(0,length(x11))))
    form = as.formula(formula$endogReg[[ii]])
    mf = model.frame(formula = form)
    m = dim(mf)[2]
    x = model.matrix(attr(mf, "terms"), data=mf)
    x1temp = as.matrix(x[,1:(m-j)]); ztemp = as.matrix(x[,(m-j+1):m])
    x2[,ii] = t(x1temp%*%pi1[[ii]] + ztemp%*%pi2[[ii]] + v[,ii])
    assign(nams[ii], x2[,ii])
  }
  colnames(x2) = nams

  # then get the means for the y regressions
  y1 = c(rep(0,n))
  form = formula$formula
  mf = model.frame(formula = as.formula(form))
  m = dim(mf)[2]
  x = model.matrix(attr(mf,"terms"),data = mf)
  x1Y = as.matrix(x[,1:(m-k)]); x2Y = as.matrix(x[,(m-k+1):m])
  y0star = x1Y%*%as.matrix(gamma1) + x2Y%*%as.matrix(gamma2) + eta
  
  
  if(family == 'lognormal'){
    if(clump == T){
      cl = sample(c(0,1),n,replace=T,prob = c(.8,.2))
      logy1star = cl*.5 + (1-cl)*(x1Y%*%as.matrix(beta1) + x2Y%*%as.matrix(beta2) +u)
      #logy1star = cl*.5 + (1-cl)*(t(beta1)%*%x1 + x2*beta2+ u) 
    }
    
    else{
      #logy1star = t(beta1)%*%x1 + t(beta2)%*%x2+ u
      logy1star = x1Y%*%as.matrix(beta1) + x2Y%*%as.matrix(beta2) + u
    }
    y1t = exp(logy1star)
  }
  else if(family == 'cragg1'){
    require(truncnorm)
    #y1star = t(beta1)%*%x1 + t(beta2)%*%x2+ u
    y1star = x1Y%*%beta1 + x2Y%*%beta2 + u
    y1t = t(as.matrix(rtruncnorm(n,a=0,mean = x1Y%*%as.matrix(beta1) + x2Y%*%as.matrix(beta2), sd = sig_u),nrow=1))
  }
  else{
    print('family must be lognormal or cragg1')
    return(None)
  }
  
  #Censor y1
  y1 = c(as.numeric(y0star>0)*y1t)
  dat = data.frame(x1, x2, z, y1 = y1)
  
  if(rho!=F){true = list(sig_u = sig_u, sig_v = sig_v, 
                         tau0 = tau0, tau1=tau1,rho=rho, 
                         beta1 = beta1,beta2 = beta2, 
                         gamma1 = gamma1,gamma2 = gamma2,
                         pi1=pi1,pi2 = pi2)}
  else{true = list(sig_u = sig_u, sig_v = sig_v, 
                   tau0 = tau0, tau1=tau1,
                   beta1 = beta1,beta2 = beta2, 
                   gamma1 = gamma1,gamma2 = gamma2,
                   pi1=pi1,pi2 = pi2)}
  adds = list( Sig = Sig, Sig_err = Sig_err, A=A,
               varx2 = var(x2), vary0star = var(y0star) )
  out = list( dat = dat, parameters = true, additional = adds )
  
  return(out)
  
}