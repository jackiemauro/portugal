###############################################################
##
## A function which outputs the covariance matrix
# 
## set y1 to the log of y1 if you want to run a loglinear model
## if you leave endogReg empty, will default to using all exog & inst
## in the first stage regressions

#issue: the input will not be a list
###############################################################

hurdleIV.covMat <- function(params){
  # make sig_err
  a = rbind(params$tau0,params$tau1)
  b = matrix(c(1,params$rho,params$rho,params$sig_u),ncol = 2, byrow = T)
  c = matrix(diag(params$j)*(unlist(params$sig_v))^2, nrow = 3)
  
  Sig_err = rbind(cbind(b,a),cbind(t(a),c))
  
  if(min(eigen(Sig_err)$values)<=0){return(Inf)}
  if((params$sig_u<=0)|(any(params$sig_v<=0))){return(Inf)}
  
  #Transformation matrix from (eta, u, v) to (y0*, log y1*, x2) (without the mean)
  A = diag(params$j+2)
  A[1,] = c(1,0,params$beta2)
  A[2,] = c(0,1,params$gamma2)
  
  Sig = A%*%Sig_err%*%t(A)
  
  Sig
  }