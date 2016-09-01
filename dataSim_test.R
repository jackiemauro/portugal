#################################################################
# Combining data simulation code into one
# To do: simulate according to a supplied formula
#         namely, need this to simulate x2's not necessarily
#         using all of the x1's
#################################################################


gen_hurdleSim <- function(n=10000,
                          family,
                          rho=F,
                          params,
                          het = FALSE,
                          clump = FALSE,
                          formula){
  
  if(family == 'cragg1'){require(MASS)}
  
  mux1 = params$mux1; sigx1 = params$sigx1
  muz = params$muz; sigz = params$sigz
  beta1 = params$beta1; gamma1 = params$gamma1; pi1 = params$pi1
  beta2 = params$beta2; gamma2 = params$gamma2; pi2 = params$pi2
  sig_v = params$sig_v; sig_u = params$sig_u
  tau0 = params$tau0; tau1 = params$tau1
  m = length(mux1) +1
  j = length(muz)
  k = length(sig_v)
  
  if(het == TRUE){
    # FIGURE OUT HOW TO PROGRAM IN HETEROSCEDASTICITY
  }
  
  Sig_x2 = matrix(diag(j)*(unlist(sig_v))^2, nrow = j)
  
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
  #Generate a matrix of the eta,u,v errors and split them out
  euv = mvrnorm(n=n, mu = c(rep(0,dim(Sig_err)[1])), Sigma = Sig_err)
  eta = euv[,1]
  u = euv[,2]
  v = matrix(c(euv[,3:dim(euv)[2]]),ncol = k)
  
  #Construct the variables
  #x1 <<- matrix(c(rep(NA,m*n)),nrow=m)
  x11<<-c(rep(1,n)) #intercept term
  #x1[1,] <<- x11
  for(ii in 2:m){
    nam <- paste("x1", ii, sep = "")
    assign(nam, c(rnorm(n,mux1[(ii-1)],sigx1[(ii-1)])), envir = .GlobalEnv)
    #x1[ii,] <<- eval(nam)
  }
  
  #z <<- matrix(c(rep(NA,j*n)),nrow=j)
  for(ii in 1:j){
    nam <- paste("z", ii, sep = "")
    assign(nam, c(rnorm(n,muz[ii],sigz[ii])), envir = .GlobalEnv)
    #z[ii,] <<- c(rnorm(n,muz[ii],sigz[ii]))
  }
  
#   # this needs to depend on the model input
#   x2 <<- matrix(c(rep(NA,j*n)),nrow=j)
#   for(ii in 1:j){
#     x2[ii,] <<- t(pi1)%*%x1 + t(pi2)%*%z + v[,ii]
#   }
    
  # problem: this needs to take the right pi's. needs to split them somehow.
  x2 <<- matrix(c(rep(NA,j*n)),nrow=j)
  for(i in 1:j){
    nam <- paste("x2", i, sep = "")
    assign(nam, c(rep(0,length(x11))), envir = .GlobalEnv)
    formula = regs$endogReg[[i]]
    mf = model.frame(formula = formula)
    m = dim(mf)[2]
    x <- model.matrix(attr(mf, "terms"), data=mf)
    x1temp = x[,1:(m-j)]; ztemp = x[,(m-j+1):m]
    x2[ii,] = t(x1temp%*%pi1[[i]] + ztemp%*%pi2[[i]])
  }
 #y0star = t(gamma1)%*%x1 + t(gamma2)%*%x2 + eta

  # then get the means for the y regressions
  y1 <<- c(rep(0,n))
  formula = regs$formula
  mf = model.frame(formula = formula)
  m = dim(mf)[2]
  x <- model.matrix(attr(mf,"terms"),data = mf)
  x1 = x[,1:(m-j)]; x2 = x[,(m-j):(m-1)]
  y0star = x1%*%gamma1 + x2%*%gamma2 + eta
  
  
  if(family == 'lognormal'){
    if(clump == T){
      cl = sample(c(0,1),n,replace=T,prob = c(.8,.2))
      logy1star = cl*.5 + (1-cl)*(x1%*%beta1 + x2%*%beta2 +u)
      #logy1star = cl*.5 + (1-cl)*(t(beta1)%*%x1 + x2*beta2+ u) 
    }
    
    else{
      #logy1star = t(beta1)%*%x1 + t(beta2)%*%x2+ u
      logy1star = x1%*%beta1 + x2%*%beta2 + u
    }
    y1t = exp(logy1star)
  }
  else if(family == 'cragg1'){
    require(truncnorm)
    #y1star = t(beta1)%*%x1 + t(beta2)%*%x2+ u
    y1star = x1%*%beta1 + x2%*%beta2 + u
    y1t = t(as.matrix(rtruncnorm(n,a=0,mean = x1%*%beta1 + x2%*%beta2, sd = sig_u),nrow=1))
  }
  else{
    print('family must be lognormal or cragg1')
    return(None)
  }
  
  #Censor y1
  y1 <<- as.numeric(y0star>0)*y1t
  
  
  if(rho!=F){true <<- list(sig_v = sig_v, sig_u = sig_u, rho=rho, tau0 = tau0, tau1=tau1,
                           beta1 = beta1,beta2 = beta2, gamma1 = gamma1,gamma2 = gamma2,
                           pi1=pi1,pi2 = pi2)}
  else{true <<- list(sig_v = sig_v, sig_u = sig_u, tau0 = tau0, tau1=tau1,
                     beta1 = beta1,beta2 = beta2, gamma1 = gamma1,gamma2 = gamma2,
                     pi1=pi1,pi2 = pi2)}
  
}