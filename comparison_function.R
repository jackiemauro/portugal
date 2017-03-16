############################################################################
#                                                                          #
#                  function comparing models for paper                     #
#                                                                          #
############################################################################

# for the paper you need to run:
# compare_mods(params = tob_pars,type="Tobit",tobit_assmp = T)
# compare_mods(params = cr_pars,type="Cragg1")
# compare_mods(params = cr_pars,type="Cragg1", inst_weak = T)

compare_mods <- function(params, type, tobit_assmp = F, inst_weak = F, true_start = F){
  source("data_simulations.R")
  source("start_values.R")
  source("loglik_lgnorm.R")
  source("loglik_cragg.R")
  source("loglik_ivtobit.R")
  require(MASS)
  require(truncnorm)
  require(xtable)
  root = "C:/Users/jackie/Dropbox/portugal/paper/draft_2016-04/JM methods paper/"
  nm.rt = type
  cap.rt = type
  
  if(tobit_assmp == T){
    # to get tobit assumptions to hold, you need:
    # sigu =1, tau1=0, 
    # beta1 = gamma1/(1-tau0/sigv^2), 
    # beta2 = (gamma2 - tau0/sigv^2) / (1-tau0/sigv^2)
    nm.rt = paste(nm.rt,"_tobAssmp",sep="")
    cap.rt = paste(cap.rt, ": Tobit Assumptions Hold",sep = "")
    
    params$tau1 <- 0
    params$sig_u <- 1
    params$beta1 <- params$gamma1/(1 - params$tau0/params$sig_v^2)
    params$beta2 <- (params$gamma2- params$tau0/params$sig_v^2)/(1 - params$tau0/params$sig_v)
  }
  
  if(inst_weak == T){
    # to weaken the instrument, have it weakly correlated with x2
    # and have there be considerable other variance in x2
    
    params$pi2 <- params$pi2/4
    params$sig_v <- params$sig_v*3
    
    nm.rt = paste(nm.rt,"_weakInst",sep="")
    cap.rt = paste(cap.rt,": Weak Instrument", sep = "")
  }
  
  # generate the dataset
  if(type == "Lognormal"){gen_lgiv(n=10000, params = params)}
  else if(type == "Tobit"){gen_tob(n=10000, params = params)}
  else if(type == "Cragg1"){gen_cragg1(n=10000, params = params)}
  else{
    print("type options are : Lognormal, Tobit and Cragg1")
    stop
  }
  
  ## Solve with all three models
  # use true start values
  if(true_start ==T){
    obs_tob = list(sig_v = params$sig_v, sig_u = params$sig_u,
               tau1 = params$tau1, beta1 = params$beta1,beta2 = params$beta2,
               pi1 = params$pi1, pi2 = params$pi2)
    obs_cr = list(sig_v = params$sig_v, sig_u = params$sig_u,
                   tau0 = params$tau0, tau1 = params$tau1, 
                   beta1 = params$beta1,beta2 = params$beta2,
                   gamma1 = params$gamma1, gamma2 = params$gamma2,
                   pi1 = params$pi1, pi2 = params$pi2)
  }
  # or get start values from naive regressions
  else{
    m = length(params$mux1)
    start_values(x1=x1,x2=x2,z=z,y1=y1,m=m)
  }
  
  # run likelihood functions
  solved_tob = optim(unlist(obs_tob), loglik_tobitiv_norho, hessian=TRUE, control = list(maxit = 5000))
  solved_cr1 = optim(unlist(obs_cr), loglik_cragg_mat, hessian=TRUE, control = list(maxit = 5000))
  
  # cragg
  cry0_pred = as.numeric((solved_cr1$par["gamma1"]*x1[1,] + solved_cr1$par["gamma2"]*x2) >0)
  cry1_star = solved_cr1$par["beta1"]*x1[1,] +solved_cr1$par["beta2"]*x2 
  cry1_pred = cry1_star[cry1_star>0]*cry0_pred[cry1_star>0]
  
  # tobit
  toby1_pred = (solved_tob$par["beta1"]*x1[1,]+
                  solved_tob$par["beta2"]*x2)
  toby1_pred = as.numeric(toby1_pred>0)*toby1_pred
  
  ### Output fit figures
  name = eval(paste(root,"Figures/",nm.rt,"_compare.png",sep=""))
  
  png(filename = name, type = 'cairo')
  par(mfrow = c(1,2), mar = c(4,3,3,1), cex = 1)
  print({
    plot(cry1_pred~y1[which(cry1_star>0)], main = "Cragg1 IV fit", xlab = "Truth",ylab = "Predicted",
         xlim = c(0,max(max(c(y1,cry1_pred)))),
         ylim = c(0,max(max(c(y1,cry1_pred)))))
    abline(0,1,col="red")
    
    plot(c(toby1_pred)~c(y1), main = "Tobit IV fit", xlab = "Truth", ylab = "Predicted",
         xlim = c(0,max(max(c(y1,toby1_pred)))),
         ylim = c(0,max(max(c(y1,toby1_pred)))))
    abline(0,1,col="red")
  })
  dev.off()
  
  par(par.old)
  
  
  ## Output tables of zeros and results
  Zname = eval(paste(root,nm.rt,"_zeroes.tex",sep=""))
  Zcap = eval(paste(cap.rt,"--Zero Detection",sep=""))
  Rname = eval(paste(root,nm.rt,"_results.tex",sep=""))
  Rcap = eval(paste(cap.rt,"--Results",sep=""))
  
  ## Compare how well they did getting zeroes right
  zeroes <- rbind("true zeroes detected" = cbind(length(which(cry0_pred>0 & y1>0))/length(which(y1>0)),
                                                 length(which(toby1_pred>0 & y1>0))/length(which(y1>0))),
                  ">0 detected" = cbind(length(which(cry0_pred==0 & y1==0))/length(which(y1==0)),
                                        length(which(toby1_pred==0 & y1==0))/length(which(y1==0)))
  )
  colnames(zeroes) = c("Cragg", "Tobit")
  rownames(zeroes) = c("True zeroes detected", "True >0 detected")
  z_table <- xtable(zeroes,caption = Zcap,label = Zname)

  
  print.xtable(z_table, type="latex", file=Zname) 

  ## get standard errors
  err_tob = sqrt(diag(solve(solved_tob$hessian)))
  err_cr = sqrt(diag(solve(solved_cr1$hessian)))

  res.tab <- matrix(c(round(params$beta1,2), round(params$beta2,2),round(params$gamma1,2),round(params$gamma2,2),
                      round(solved_tob$par["beta1"],2),round(solved_tob$par["beta2"],2),NA,NA,
                      round(err_tob["beta1"],2), round(err_tob["beta2"],2),NA,NA,
                      round(solved_cr1$par["beta1"],2),round(solved_cr1$par["beta2"],2),round(solved_cr1$par["gamma1"],2),round(solved_cr1$par["gamma2"],2),
                      round(err_cr["beta1"],2), round(err_cr["beta2"],2),round(err_cr["gamma1"],2),round(err_cr["gamma2"],2)
                      ),
                    ncol = 4, byrow = T)
  colnames(res.tab)<-c("Beta1", "Beta2","Gamma1","Gamma2")
  rownames(res.tab)<-c("True Values", "Tobit IV", "Tobit SD","Cragg 1 IV", "Cragg SD")

  r_table <- xtable(res.tab,caption = Rcap, label = Rname)
  res.xtab <- print.xtable(r_table, type = "latex", file = Rname)
  

  ## Final output
  comp.out <<- list(tobit = solved_tob,
                    cragg = solved_cr1,
                    true = unlist(params))
  
}


