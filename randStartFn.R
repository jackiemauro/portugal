# run at random start points
rand.Start<- function(n, silent = F){
  start = hurdleIV.start_vals(regs, family = 'lognormal')
  start = start[-5] # get rid of rho
  start = tagBeg(start)
  randOut = matrix(c(rep(NA,n*length(start))),ncol = length(start))
  converge = NA; ll = NA
  randEig =  matrix(c(rep(NA,n*length(start))),ncol = length(start))
  randStart = matrix(c(rep(NA,n*length(start))),ncol = length(start))
  
  for(i in 1:n){
    if(silent == F){print(paste("round ",i," of ",n, sep = ""))}
    val = sapply(start, function(x) x + rnorm(1,0,sd = abs(x*2)))
    val['sig_u_ls1.elem1'] = abs(val['sig_u_ls1.elem1'])
    val['sig_v_ls2.elem1.subelem1'] = abs(val['sig_v_ls2.elem1.subelem1'])
    err = 0
    tryCatch({
      result = optim(val,loglik_lgnorm
            , control = list(maxit=5000)
            , hessian = T
      )
      randStart[i,] = val
      randOut[i,] <- result$par
      randEig[i,] <- eigen(result$hessian)$value
      converge[i] <- result$convergence
      ll[i] <- result$value
      if(min(eigen(result$hessian)$value)<0 & silent == F){
        print(paste("BAD HESSIAN: round",i))
      }
    }, error = function(e){
      print(paste("optimizer error, skipping round ",i,sep = ""))
      randStart[i,] = val
      randOut[i,] <- rep(NA,length(start))
      randEig[i,] <- rep(NA,length(start))
      converge[i] <- NA
      ll[i] <- NA
    }
    )
  }
  all.out = list(startvals = randStart
                 , estimates = randOut
                 , eigenvalues = randEig
                 , convergence = converge
                 , loglikelihood = ll)
  return(all.out)
}
