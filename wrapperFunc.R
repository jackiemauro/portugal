hurdle.IV<-function(formula,
                    inst,
                    endog,
                    exog,
                    data,
                    endog_reg = F,
                    start_val = F,
                    options = list(cholesky = T
                                   , maxit = 5000
                                   , trace = 0
                                   , method = "BFGS")
                    ){
  
  ###### format your data, grouping variable types #####
  attach(data)
  inst_mat = as.data.frame(matrix(inst,nrow = dim(data)[1], byrow = F))
  names(inst_mat)<-rename.input(substitute(inst))
  
  endog_mat = as.data.frame(matrix(endog, nrow = dim(data)[1], byrow = F))
  names(endog_mat)<-rename.input(substitute(endog))
  
  exog_mat = as.data.frame(matrix(exog, nrow = dim(data)[1], byrow = F))
  names(exog_mat)<-rename.input(substitute(exog))
  
  outcome = eval(parse(text=paste("data$",formula[[2]],sep = "")))
  y_mat = model.matrix(formula)
  y_colnames = colnames(y_mat)
  y_exog = exog_mat[names(exog_mat) %in% y_colnames]
  y_endog = endog_mat[names(endog_mat) %in% y_colnames]
  if(any(names(inst_mat)%in% y_colnames)){
    print("Error: main regression cannot include instruments")
    return(NA)
  }
  
  
  ######### checks #########
  #check you have a formula
  tryCatch({
    form = as.formula(formula)
    }, error = function(e){
      print("Error: formula must be a formula, eg: y~x1+x2")
      return(NA)
    }
    )
  # check the formula includes endogenous variables
  if(dim(y_endog)[2]==0){
    print("Error: endogenous variables must be included in main regression")
    return(NA)
  }
  
  # check no more endogenous variables than instruments
  if(length(inst)<length(endog)){
    print("Error: More endogenous variables than instruments")
    return(NA)
  }
  
  ##### make endog reg #####
  # if not specified, replace endog_reg with formula that has all exog's, all inst's
  if(length(endog_reg) == 0){
    a = names(endog_mat)[1]
    b = paste(names(exog_mat),collapse = "+")
    d = paste(names(inst_mat), collapse = "+")
    endog_reg = list(as.formula(paste(a,"~",b,"+",d)))
  }
  else{
    if(length(endog_reg)!=dim(endog_mat)[2]){
      print("Error: Need one regression for each endogenous variable")
      return(NA)
    }
  }
  
  ER_mat = lapply(endog_reg, function(x) model.matrix(x))
  ER_colnames = lapply(ER_mat, function(x) colnames(x))
  ER_exog = lapply(ER_colnames, function(x) exog_mat[names(exog_mat) %in% x])
  ER_inst = lapply(ER_colnames, function(x) inst_mat[names(inst_mat) %in% x])
  
  ### drop NA's
  mf = cbind(outcome, exog_mat, inst_mat, endog_mat)
  mf = mf[complete.cases(mf),] #drop NA's
  detach(data)
  
  # with new datasets, start from here
  ############# get start values #######
  # if start values aren't specified, get start values
  if(start_val == F){
    start_val = start.val(formula = update(form,outcome~.)
                          , endog_reg = endog_reg
                          , data = mf)
  }
  
  else{
    if(class(start_val) != "list"){
      print("Error: start values must be a list")
      return(NA)
    }
    if(length(start_val[['beta']]!=length(start_val[['gamma']]))){
      print("Error: beta and gamma must be equal length vectors")
      return(NA)
    }
    k = dim(endog_mat)[2]
    if(length(start_val$pi)!=k | length(start_val$tau0)!= k | length(start_val$tau1)!=k){
      print("Error: tau0, tau1 and number of pi coordinates must equal number of endogenous variables")
      return(NA)
    }
  }
  
  # start cov values should be cholesky-fied if that option is true
  cov_start = make.cov(rho=start_val$rho
                       ,tau0=start_val$tau0
                       ,tau1=start_val$tau1
                       ,y_sd=start_val$y_sd
                       ,endog_sd=start_val$endog_sd)
  if(options$cholesky == T){
    cov_start = chol(cov_start)
  }

  
  ########## run optimizer #############
  out = optim(start_vec
              , loglik_lgiv
              , method= options$method
              , hessian  = T
              , control = list(maxit = options$maxit,trace = options$trace)
  )

}

rename.input <- function(input){
  if(class(input)=="name"){
    name = toString(input)
  }
  else if(class(input) == "call"){
    text = input[-1]
    k = length(text)
    name = NULL
    for(i in 1:k){name[i] = toString(text[[i]])}
    return(name)
  }
}

make.cov<- function(rho,tau0,tau1,y_sd,endog_sd){
      mat1 = matrix(c(1,rho,rho,y_sd^2),ncol = 2, byrow = T)
      tau_mat = as.matrix(cbind(unlist(tau0),unlist(tau1)))
      endog_mat = diag(length(endog_sd))*endog_sd^2
      Sig_err = rbind(cbind(mat1,t(tau_mat)),cbind(tau_mat,endog_mat))
      
      return(as.matrix(Sig_err))
    }



