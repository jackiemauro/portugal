# use usual iv regressions and bootstrap errors
require(boot)

sep_est <- function(type,y,endog,exog,inst,data){
  d<-data[indices,]
  require('AER')
  require('ivprobit')
  attach(data)
  
  y0 = as.numeric(y>0)
  
  lin = lm(endog~exog + inst)
  prob = ivprob(y=y0,x1=exog,y2=endog,x=data.frame(exog,inst))
  if(type == "cragg"){
    reg = ivreg(y ~ exog + endog | exog + inst, subset = y>0)
  }
  else{
    reg = ivreg(log(y) ~ exog + endog | exog + inst, subset = y>0)
  }
  
  
  coefs = data.frame(firststage = coef(lin)
                     ,probit =  coef(prob)
                     ,secondstage = coef(reg))

  detach(data)
  return(coefs)
}

sep_est_boot <- function(data,indices){
  d<-data[indices,]
  require('AER')
  require('ivprobit')
  
  y0 = as.numeric(y>0)
  
  lin = lm(endog~exog1 + inst1, data = d)
  prob = ivprob(y=d$y0,x1=d$exog1,y2=d$endog,x=data.frame(d$exog1,d$inst1))
  if(type == "cragg"){
    reg = ivreg(y ~ exog1 + endog | exog1 + inst1,data = d, subset = y>0)
  }
  else{
    reg = ivreg(log(y) ~ exog1 + endog | exog1 + inst1,data = d, subset = y>0)
  }
  
  
  #   coefs = data.frame(firststage = coef(lin)
  #                      ,probit =  coef(prob)
  #                      ,secondstage = coef(reg))
  
  coefs = c(coef(lin),c(coef(prob)), coef(reg))
  
  return(coefs)
}

type = "cragg"
boot.results <- boot(data=dat,statistic=sep_est_boot,R=100)

