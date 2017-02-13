hurdle.IV<-function(formula,
                    inst,
                    endog,
                    exog,
                    data,
                    endog_reg = F,
                    start_val = list(),
                    type = "lognormal",
                    options = list(cholesky = T
                                   , maxit = 5000
                                   , trace = 0
                                   , method = "BFGS")
                    ){
  
  ###### format your data, grouping variable types #####
  attach(data)
  inst_mat = as.data.frame(matrix(inst,nrow = dim(data)[1], byrow = F))
  inst_names = rename.input(substitute(inst))
  names(inst_mat)<-inst_names
  
  endog_mat = as.data.frame(matrix(endog, nrow = dim(data)[1], byrow = F))
  endog_names = rename.input(substitute(endog))
  names(endog_mat)<- endog_names
  
  exog_mat = as.data.frame(matrix(exog, nrow = dim(data)[1], byrow = F))
  exog_names = rename.input(substitute(exog))
  names(exog_mat)<- exog_names
  
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
      detach(data)
      return(NA)
    }
    )
  # check the formula includes endogenous variables
  if(dim(y_endog)[2]==0){
    print("Error: endogenous variables must be included in main regression")
    detach(data)
    return(NA)
  }
  
  # check no more endogenous variables than instruments
  if(length(inst)<length(endog)){
    print("Error: More endogenous variables than instruments")
    detach(data)
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
      detach(data)
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
  attach(mf)
  
  ############# get start values #######
  # if start values aren't specified, get start values
  if(length(start_val) == 0){
    start_val = start.val(formula = update(form,outcome~.)
                          , endog_reg = endog_reg
                          , data = mf
                          , type = type)
  }
  
  else{
    if(class(start_val) != "list"){
      print("Error: start values must be a list")
      detach(mf)
      return(NA)
    }
    if(length(start_val[['beta']])!=length(start_val[['gamma']])){
      print("Error: beta and gamma must be equal length vectors")
      detach(mf)
      return(NA)
    }
    k = dim(endog_mat)[2]
    if(length(start_val$pi)!=k | length(start_val$tau0)!= k | length(start_val$tau1)!=k){
      print("Error: tau0, tau1 and number of pi coordinates must equal number of endogenous variables")
      detach(mf)
      return(NA)
    }
    if(is.null(names(start_val$gammas))){
      names(start_val$gamma)<-c("Intercept", exog_names, endog_names)
#       num_gam = length(start_val$gamma)
#       names_gam = c("Intercept", strsplit(toString(formula[[3]]),", ")[[1]][2:num_gam])
#       names(start_val$gamma)<-names_gam
    }
    if(is.null(names(start_val$beta))){
      names(start_val$beta)<-c("Intercept", exog_names, endog_names)
#       num_bet = length(start_val$beta)
#       names_bet = c("Intercept", strsplit(toString(formula[[3]]),", ")[[1]][2:num_bet])
#       names(start_val$beta)<-names_bet
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
  # save info about parameters and model matrices
  pars = list(len_cov = length(which(upper.tri(cov_start, diag = T)))-1
              ,num_endog = dim(endog_mat)[2]
              ,num_betas = dim(y_mat)[2]
              ,num_pis = lapply(ER_mat, function(x) dim(x)[2])
              ,myChol = options$cholesky
              ,ER_mat = ER_mat
              ,y_mat = y_mat
              ,endog_mat = endog_mat
              ,inst_names = inst_names
              ,exog_names = exog_names
              ,endog_names = endog_names)
  attach(pars)
  
  cov_in = cov_start[upper.tri(cov_start, diag = T)][-1]
  start_vec = c(cov_in, start_val$beta, start_val$gamma, unlist(start_val$pi))
  
  if(type == "lognormal"){
    ll = loglik_lgiv
  }
  else if(type == "cragg"){
    ll = loglik_craggiv
  }

  out = optim(start_vec
              , ll
              , method= options$method
              , hessian  = T
              , control = list(maxit = options$maxit,trace = options$trace)
  )
  name_pars = name.pieces(out$par)
  name_pars$cov = reconstitute.cov(vals=name_pars$cov_start
                                   ,num=num_endog
                                   ,chol=myChol)
  
  detach(mf)
  detach(pars)
  
  if(min(eigen(out$hessian)$value)*max(eigen(out$hessian)$value)<0){
    print("bad hessian")
  }
  return(list(out,name_pars))
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


make.covTrans <- function(a,num_endog,gamma,beta,option = "mat",noname = F){
  #option "parameters": a = list(rho,tau0,tau1,y_sd,endog_sd)
  #option "mat": a is full matrix
  #option "vector": a is vector of all elements of matrix
  
  if(option == "mat"){
    if(is.matrix(a)){
      Sig_err = a
    }
    else{
      cat("Error: You have set option to 'mat', please supply a matrix\n
          Other options:\n
          option 'vector': Insert full matrix as a vector\n
          option 'parameters': Insert only the parameters of interest as list, unimportant off-diagonal elements will be set to 0")
      return(NA)
    }
  }
  
  else if(option=="vector"){
    tryCatch({
      Sig_err = matrix(a, ncol = 2+num_endog, byrow = F)
    }, error = function(e){
      cat("Error: You have not supplied the correct length vector \n
          Vector should be (2+number of endogenous variables)^2 long\n
          Other options:\n
          option 'mat': Insert full matrix\n
          option 'parameters': Insert only the parameters of interest as list, unimportant off-diagonal elements will be set to 0")
      return(NA)
    }, warning = function(w){
      cat("Error: You have not supplied the correct length vector \n
          Vector should be (2+number of endogenous variables)^2 long\n
          Other options:\n
          option 'mat': Insert full matrix\n
          option 'parameters': Insert only the parameters of interest as list, unimportant off-diagonal elements will be set to 0")
      return(NA)
    }
      )
    
    }
  
  else if(option == "parameters"){
    tryCatch({
      mat1 = matrix(c(1,a$rho,a$rho,a$y_sd^2),ncol = 2, byrow = T)
      tau_mat = as.matrix(cbind(a$tau0,a$tau1))
      endog_mat = diag(length(a$endog_sd))*a$endog_sd^2
      Sig_err = rbind(cbind(mat1,t(tau_mat)),cbind(tau_mat,endog_mat))
    }, error = function(e){
      print(e)
      cat("Error: You have not supplied the correct length vector \n
          Vector should include rho, y_sd, tau0, tau1, endog_sd\n
          Other options:\n
          option 'mat': Insert full matrix\n
          option 'vector': Insert full matrix as a vector\n")
      return(NA)
    }, warning = function(w){
      print(w)
      cat("Error: You have not supplied the correct length vector \n
          Vector should include rho, y_sd, tau0, tau1, endog_sd\n
          Other options:\n
          option 'mat': Insert full matrix\n
          option 'vector': Insert full matrix as a vector\n")
      return(NA)
    }
      )
    }
  
  
  A = diag(2+num_endog)
  
  if(noname == T){
    gam2 = tail(gamma,num_endog)
    bet2 = tail(beta,num_endog)
  }
  else{
    gam2 = gamma[names(gamma) %in% endog_names]
    bet2 = beta[names(beta) %in% endog_names]
  }
  
  A[1,3:(2+num_endog)] <- gam2
  A[2,3:(2+num_endog)] <- bet2
  
  Sig = A%*%Sig_err%*%t(A)
  
  if(min(eigen(Sig)$values)<=0){
    print("bad start sigma values")
    return(NA)
  }
  
  return(Sig)
  }

reconstitute.cov<-function(vals,num,chol=myChol){
  # de-cholesky-ify if you're supposed to
  if(chol == T){
    cov_vals = c(1,vals)
    empty = diag(2+num)
    empty[upper.tri(empty,diag = T)]<-cov_vals
    Sig_err = t(empty)%*%empty
    Sig_err = Sig_err/Sig_err[1,1]
  }
  else{
    print("Warning: no cholesky decomposition makes maximization less stable")
    cov_vals = c(1,vals)
    empty = diag(2+num)
    empty[upper.tri(empty,diag = T)]<-cov_vals
    tempty = t(empty)
    tempty[upper.tri(tempty,diag = T)]<-cov_vals
    Sig_err = tempty
  }
  return(Sig_err)
}

cant.solve <- function(m) class(try(solve(m),silent=T))!="matrix"

name.pieces<-function(t){  
  cov_start = t[1:len_cov]
  beta = t[(len_cov+1):(len_cov+num_betas)]
  gamma = t[(len_cov+num_betas+1):(len_cov+2*num_betas)]
  pis = list()
  k = len_cov+2*num_betas+1
  for(i in 1:num_endog){
    pis[[i]] = t[k:(k+num_pis[[i]]-1)]
    k = k+num_pis[[i]]
  }
  #pis = t[(len_cov+2*num_betas+1):length(t)]
  
  return(list(cov_start = cov_start
              ,beta = beta
              ,gamma = gamma
              ,pi = pis))
}
