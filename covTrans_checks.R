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

#comparison
t1=make.covTrans(list(rho = mat[1,2], y_sd = mat[2,2],
                   endog_sd = c(mat[3,3],mat[4,4]),
                   tau0 = c(mat[1,3:4]),tau1 = c(mat[2,3:4]))
              , num_endog = 2, gamma = c(1,2), beta = c(3,3),
              option = "parameters")

t2=make.covTrans1(rho = mat[1,2], y_sd = mat[2,2],
               endog_sd = c(mat[3,3],mat[4,4]),
               tau0 = c(mat[1,3:4]),tau1 = c(mat[2,3:4]),
               gamma = c(1,2), beta = c(3,3))

t1==t2

t1=make.covTrans(list(rho = 3, y_sd = 2,
                      endog_sd = c(2,3),
                      tau0 = c(3,4),tau1 = c(5,6))
                 , num_endog = 2, gamma = c(1,2), beta = c(3,3),
                 option = "parameters")

t2=make.covTrans1(rho = 3, y_sd = 2,
                  endog_sd = c(2,3),
                  tau0 = c(3,4),tau1 = c(5,6),
                  gamma = c(1,2), beta = c(3,3))

t1==t2v

########## checking reconstituting fn  ##############
myChol = T
mat = Posdef(4)
st = mat/mat[1,1]
c = chol(st)
t = c[upper.tri(c,diag = T)][-1]

st2 = reconstitute.cov(vals =t, num = 1)
st2 = reconstitute.cov(vals =t, num = 2)
st2 == st

#no cholesky
myChol = F
mat = matrix(c(1,2,3,
               2,4,5,
               3,5,9),
             ncol = 3, byrow = T)
t = mat[upper.tri(mat,diag = T)][-1]

st2 = reconstitute.cov(vals =t, num = 1)
st2 == mat