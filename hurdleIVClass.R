###############################################################
##
## The purpose of this file is to start developing the 
## class of objects which will be hurdle IV objects for analysis
##
###############################################################

# make it its own thing
hurdleIV <- function(x,...) UseMethod("hurdleIV")

getLik <- function(family){
  if(family == 'cragg1'){
    source("loglik_cragg.R")
    loglik = loglik_cragg
  }
  else if(family == 'lognormal'){
    source("loglik_lgnorm.R")
    loglik = loglik_lgnorm
  }
  else{
    print("Choose family from 'cragg1' or 'lognormal")
    loglik = NA 
  }
  return(loglik)
}

tagBeg <- function(b){
  # this function just tags the start of lists by 
  # appending a tag to their name
  
  for( ls in 1:length(b) ){
    names(b)[ls] <- paste(names(b)[ls],paste0('ls',ls),sep = "_")
    
    for( elm in 1:length(b[[ls]]) ){
      names(b[[ls]])[elm] <- paste0('elem',elm)
      
      for( subelm in 1:1:length(b[[ls]][[elm]]) ){
        names(b[[ls]][[elm]])[subelm] <- paste0('subelem',subelm)
      }
    }
  }
  
  return(unlist(b))
}

hurdleIV.function <- function(formula, family = 'cragg1', data = list(),...){
  
  # get starting values
  if(start_values == F){
    starts = hurdleIV.start_vals(formula)
  }
  else{
    starts = start_values
  }
  
  j = length(endog); l = length(inst)
  
  # need to figure out how to get this dataset into the right shape
  mf <- model.frame(formula=formula$formula, data=data)
  x <- model.matrix(attr(mf, "terms"), data=mf)
  y1 <- model.response(mf)
  
  loglik = getLik(family)
  
  # get likelihood
  solved = optim(unlist(start_values), loglik, hessian=TRUE)
  solved
}