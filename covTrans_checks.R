source("loglikFunc.R")

make.covTrans(5)

make.covTrans(matrix(c(1,2,3,4),ncol = 2))

make.covTrans(matrix(c(1,2,3,4,5,6,8,9),ncol = 3),num_endog = 1)

make.covTrans(matrix(c(1,2,3,4,5,6,7,8,9),ncol = 3)
              ,num_endog = 1
              ,gamma = 0
              ,beta = 0)

mat = Posdef(4)
make.covTrans(mat,num_endog = 2, gamma=0,beta=0)

# with vector options
make.covTrans(list(rho = 3, y_sd = 2, endog_sd = c(1,2), tau0 = c(3,4),
                   tau1 = c(5,6)), num_endog = 2, gamma = c(0,0), beta = c(0,0))

# with parameter options
make.covTrans(list(rho = 3, y_sd = 2, endog_sd = c(1,2),
                   tau0 = c(3,4),tau1 = c(5,6))
              , num_endog = 2, gamma = c(0,0), beta = c(0,0),
              option = "parameters")

make.covTrans(list(rho = mat[1,2], y_sd = mat[2,2],
                   endog_sd = c(mat[3,3],mat[4,4]),
                   tau0 = c(mat[1,3:4]),tau1 = c(mat[2,3:4]))
              , num_endog = 2, gamma = c(0,0), beta = c(0,0),
              option = "parameters")