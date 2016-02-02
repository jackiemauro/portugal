#################################################################################
# starting from scratch Jan 27, 2016
# file should simulate a dataset from a lognormal hurdle model
# then it should run an mle procedure to recover the parameters of the model

# currently the true values don't actually maximize this likelihood 
#################################################################################

require(stats4)
require(MASS)
require(maxLik)
require(optimx)

##################### simulate the dataset ##################

# error terms
# setting one tau high. code works ok if you set these both very small
sig_u = .5; sig_v = .7
tau1 = .02; tau0 = .5
sig = matrix( c(sig_u^2, 0,    tau1,
                0,     1,    tau0,
                tau1,  tau0, sig_v^2),
              ncol = 3, byrow = T)
uev = mvrnorm(n=n, mu = c(0,0,0), Sigma = sig)
u = uev[,1]; eta = uev[,2]; v = uev[,3]

# exogenous variables
n = 10000
x1 = rnorm(n, .5, sqrt(.02)); z = rnorm(n, .2, sqrt(.02))
Z = matrix(c(x1, z), ncol = 2, byrow = F)

# endogenous variable
Pi = c(-.2, .7)
x2 = Pi%*%t(Z) + v
X = matrix(c(x1,x2), ncol = 2, byrow = F)

# second stage regressions
Gamma = c(2, -.8); Beta = c(.02, .005)
y0s = Gamma%*%t(X) + eta
y0 = ifelse(y0s>0,1,0)

lny = Beta%*%t(X) + u
y1 = c(exp(lny)*y0)

# check look
hist(y1)
hist(log(y1[y1>0]))

############################ mle #######################

loglik <- function(sig_v, sig_u, tau1,tau0,beta1, beta2,gamma1,gamma2,pi1,pi2){
  
  ## conditional means and variances
  
  # helpers
  alpha0_2 = tau0/(sig_v^2); alpha1_2 = tau1/(sig_v^2)
  alpha0_12 = tau0/((sig_u^2)*(sig_v^2) - (tau1^2))
  
  # parameters
  mu0_12 = gamma1*x1 + gamma2*x2 + 
    alpha0_12*(-tau1*(log(y1) - (beta1*x1 + beta2*x2)) + (sig_u^2)*(x2 - (pi1*x1 + pi2*z)))
  sig0_12 = 1 - alpha0_12*tau0*sig_u^2
  if(sig0_12 <= 0){
    print("sig0_12 < 0, replacing with NA")
    sig0_12 <- NA
  } 
  
  mu0_2 = gamma1*x1 + gamma2*x2 + alpha0_2*(x2 - (pi1*x1 + pi2*z))
  sig0_2 = 1 - alpha0_2*tau0
  if(sig0_2 <= 0){
    print("sig0_2 < 0, replacing with NA")
    sig0_2 <- NA
  } 
  sig0_2 <- ifelse(sig0_2 > 0, sig0_2, NA)
  
  mu1_2 = beta1*x1 + beta2*x2 + alpha1_2*(x2 - (pi1*x1 + pi2*z))
  sig1_2 = sig_u^2 - alpha1_2*tau1 
  if(sig1_2 <= 0){
    print("sig1_2 < 0, replacing with NA")
    sig1_2 <- NA
  } 
  
  ## densities
  
  f1 = 1-pnorm(mu0_2/sig0_2)
  f1[y1 > 0] <- 1
  f1[is.infinite(f1)] <- NA
  
  f2 = pnorm(mu0_12/sig0_12)*dnorm(log(y1), mu1_2, sig1_2)
  f2[y1 == 0] <- 1
  f2[is.infinite(f2)] <- NA
  
  full = c(f1*f2)
  if(length(full[is.infinite(full) || is.nan(full) || full==0 || is.na(full)]) > 0){
    print("infinite/nan/0/na values of likelihood")
  }
  full[is.infinite(full) || is.nan(full) || full==0 || is.na(full)] <- 1e-6
  
  
  ## log likelihood
  
  ll = -sum( log(full) )
  print(ll)
  ll
  
}


# get starting values
x2.v <- c(x2)
t1 <- lm(x2.v ~ x1 + z)
I = ifelse(y1>0,1,0)
t2 <- glm(I~x1+x2.v, family = binomial(link = "probit"))
t3 <- lm(log(y1[y1>0])~x1[y1>0]+x2.v[y1>0])

obs = list(sig_v=sd(t1$res), sig_u=sd(t3$res), 
       tau1=0,tau0=0,beta1=t3$coef[[2]], beta2=t3$coef[[3]],
       gamma1=t2$coef[[2]],gamma2=t2$coef[[3]],
       pi1=t1$coef[[2]],pi2=t1$coef[[3]])


# use true values
true = list(sig_v = sig_v, sig_u = sig_u, tau1= tau1, tau0=tau0,
           beta1 = Beta[1], beta2 = Beta[2], gamma1= Gamma[1],
           gamma2 = Gamma[2], pi1=Pi[1], pi2 = Pi[2])

# run program
lgiv =  mle(loglik, start = true, method = "L-BFGS-B",
            lower = c(1e-6,1e-6,rep(-Inf,8)),
            upper = c(rep(Inf,10)))

# compare coefficients to truth
comp<-data.frame(mle.results = c(lgiv@coef), 
                 truth = c(sig_v, sig_u, tau1, tau0,Beta[1],Beta[2],
                           Gamma[1], Gamma[2], Pi[1], Pi[2]))
View(comp)
