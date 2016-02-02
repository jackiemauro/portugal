#######################################################################
# runs the maxLik algorithm on simulated dataset
# uses the bayes law trick Maria thought of
# uses better notation and generates endogeneity differently
# have not updated the ml function yet!
#######################################################################

require(stats4)
require(MASS)
require(maxLik)
require(optimx)

################### make a simulated dataset ####################

# error terms for all 3 regressions: model all 3 as jointly mv normal
sig_u = 1.5; sig_v = 1.3; tau1 = .8; tau0 = .6
sig = matrix(c(sig_u^2, 0,    tau1, 
               0,       1,    tau0, 
               tau1,    tau0, sig_v^2), 
             ncol =3, byrow = T)
uev = mvrnorm(n = 1000, mu = c(0,0,0), Sigma = sig)
u = uev[,1]; eta = uev[,2]; v = uev[,3]

# exogenous covariates: x1 exogenous; z instrument
x1 = rnorm(1000,8,2); z = rnorm(1000,2,1)

# first stage regression
pi1 = .3; pi2 = 2
x2 = pi1*x1 + pi2*z + v

# construct probit section
gam1 = .08; gam2 = .05
y0 = x1*gam1 + x2*gam2 + eta
pzero = 1-pnorm(y0)
s = rbinom(1000, 1, 1-pzero) #s = 1 if y is larger than 0

# unobserved variable
beta1 = .07; beta2 = .025
lny = x1*beta1 + x2*beta2 + u

# observed variable
y1 = exp(lny)*s

##################### ml procedure ###################

# define the function that should be maximized
# needs to get improved to be vectors not scalars
loglik <-function(sigma_v, sigma_u, Gamma1, Gamma2, Beta1, Beta2, tau1, tau0, pi1, pi2){
  I = y1==0                       
  
  # helper values
  alpha1 = tau1/(sigma_v^2); alpha0 = tau0/(sigma_v^2)
  mu1_2 = x1*Beta1 + x2*Beta2 + alpha1*(x2 - x1*pi1 - z*pi2)
  sig1_2 = sigma_u^2 - ((tau1^2)/(sigma_v^2))
  
  mu0_2 = x1*Gamma1 + x2*Gamma2 + alpha0*(x2 - x1*pi1 - z*pi2)
  sig0_2 = 1 - ((tau0^2)/(sigma_v^2))
  
  mu0_12 = x2*Gamma2 + x1*Gamma1 - 
    (tau1*tau0/(sigma_u^2*sigma_v^2 - tau1^2))*(log(y1) - (x1*Beta1 + x2*Beta2)) + 
    (sigma_u^2*tau0/(sigma_u^2*sigma_v^2 - tau1^2))*(x2 - (x1*pi1 + z*pi2))
  mu0_12[is.infinite(mu0_12)] <- NA
  mu0_12[is.nan(mu0_12)] <- NA
  sig0_12 = 1 - (sigma_u^2*tau0^2/(sigma_u^2*sigma_v^2 - tau1^2))
  
  # density
  f = ( (pnorm(0,mu0_2,sig0_2))^I )*
    ( (pnorm(0,-mu0_12,sig0_12) * dnorm(log(y1), mu1_2, sig1_2) * (1/(y1*sig1_2)))^(1-I) )
  #pnorm(0, -mu, sigma) = Phi( (0+mu)/sigma ) = Phi(mu/sigma)
  
  # density tweaks
  f[is.nan(f)] <- NA
  f[is.infinite(f)] <- NA
  f[f==0]<-NA
  
  # likelihood
  log.lik = -sum(log(f), na.rm = T)
  return(log.lik)
}


################ test with starting values of varying difficulty ################

# trivial: start them at true values
true = list(sigma_v = sd(v), sigma_u = sd(u), 
            Gamma1 = gam1, Gamma2 = gam2,
            Beta1 = beta1, Beta2 = beta2, 
            tau1 = cov(eta,v)^2/var(v), tau0 = cov(u,v)^2/var(v),
            pi1 = pi1, pi2 = pi2)

lgiv = mle(loglik, start = true, method = "L-BFGS-B", 
           lower = c(0,0,rep(-Inf,8)), upper = c(rep(Inf,10)))
View(cbind(lgiv@coef,true)) 

# little harder: inject some normal noise
foo <- function(x) x + rnorm(1,0,.04)
harder = lapply(true, foo)
lgiv2 =  mle(loglik, start = harder, method = "L-BFGS-B", 
             lower = c(0, 0, rep(-Inf,8)),
             upper = c(rep(Inf,10)))
View(cbind(lgiv2@coef,true))

# even harder: just set everything to 2
foo <- function(x) 2 
rlhard = lapply(true,foo)
lgiv3 =  mle(loglik, start = rlhard, method = "L-BFGS-B", 
             lower = c(0, 0,rep(-Inf,8)),
             upper = c(rep(Inf,10)))
View(cbind(lgiv3@coef,true)) #can't do it

# realistic: run naive regressions to get start values
start.lm = lm(x2~x1+z)
o.pi1 = start.lm$coef[2]; o.pi2 = start.lm$coef[3]; o.sigma_v = summary(start.lm)$sigma

bit = ifelse(y1==0,0,1)
start.probit = glm(bit ~ x2 + x1, family=binomial(link="probit"))
o.gamma1 = start.probit$coef[2]; o.gamma2 = start.probit$coef[3]

start.lg = lm(log(y1[y1>0]) ~ x2[y1>0] + x1[y1>0])
o.beta1 = start.lg$coef[3]; o.beta2 = start.lg$coef[2]

obs.start = list(sigma_v = as.numeric(o.sigma_v), sigma_u =1, 
                 Gamma1 = as.numeric(o.gamma1), Gamma2 =as.numeric(o.gamma2),
                 Beta1 = as.numeric(o.beta1), Beta2 = as.numeric(o.beta2), 
                 tau1 = 0, tau0 = 0, 
                 pi1 = as.numeric(o.pi1), pi2 = as.numeric(o.pi2))

mu0_2 = x2*as.numeric(o.gamma2) + x1*as.numeric(o.gamma1) 
sig0_2 = 1

lgiv.o =  mle(loglik, start = obs.start, method = "L-BFGS-B", 
              lower = c(0, 0,0,rep(-Inf,8)),
              upper = c(rep(Inf,11)))
View(cbind(lgiv.o@coef,obs.start,true)) 
