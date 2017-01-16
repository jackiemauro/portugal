require(MASS)
require(maxLik)


# make a simulated dataset that will do what we want
x1 = rnorm(1000,8,2)
x2 = rnorm(1000,2,1)

pi1 = .3
pi2 = 2

sig = matrix(c(1.1,.2,.2,.04), ncol =2, byrow = T)
uv = mvrnorm(n = 1000, mu = c(0,0), Sigma = sig)
u = uv[,1]
v = uv[,2]

y2 = x1*pi1 + x2*pi2 + v

gam1 = -.80
gam2 = .5

beta1 = .02
beta2 = .025

eta = rnorm(1000,0,1)
y0 = y2*gam1 + x1*gam2 + eta
pzero = 1-pnorm(y0)
s = rbinom(1000, 1, 1-pzero)

lny = y2*beta1 + x1*beta2 + u
y1 = exp(lny)*s

I = y1==0


# define the function that should be maximized
loglik <-function(params){
  # indicator for above the cutoff tau=0
  I = y1==0                       
  
  # starting values in maxLik
  sigma_v = params[1]; sigma_uv = params[2]; sigma_ev = params[3]
  Gamma1 = params[4]; Gamma2 = params[5];
  Delta1 = params[6]; Delta2 = params[7]
  alpha1 = params[8]; alpha2 = params[9]
  pi1 = params[10]; pi2 = params[11]
  
  # helper values
  b = y1*Delta1 + x1*Delta2 + alpha1*(y2 - x1*pi1 - x2*pi2)
  m = (y1*Gamma1 + x1*Gamma2 + alpha2*(y2 - x1*pi1 - x2*pi2))/sigma_ev
  temp = pnorm(m)*(( 1/( 2*pi*sigma_uv ) )*exp(-( 1/(2*sigma_uv) ) * (log(y1) - b)^2))
  
  # break out likelihoods
  f2 = -0.5*log(2*pi) - log(sigma_v) - 0.5*(y2 - x1*pi1 - x2*pi2)^2/sigma_v^2
  f1 = I*log( 1 - pnorm(m) ) + 
    (1-I)*log(temp)
  f1[is.na(f1)] <- 0
  f1[is.nan(f1)] <- 0
  f1[is.infinite(f1)] <- 0
  
  # log likelihood fn for iv lognormal hurdle
  ll = f1 + f2
  ll                                                
}
#has log(pnorm) values which go to negative inf if y1 is too big

alpha1 = cov(u,v)/var(v)
alpha2 = cov(eta,v)/var(v)


true = c(sd(v), sqrt(var(u) - cov(u,v)^2/var(v)), 
         sqrt(1 - cov(eta,v)^2/var(v)), gam1,gam2,beta1,beta2, 
         alpha1,alpha2,pi1,pi2)
reslg1 = maxLik(loglik, start = true)

harder = true + rnorm(1,0,.04)
reslg2 = maxLik(loglik, start = harder)

names = c("sv","suv","sev","g1","g2","b1","b2","a1","a2","p1","p2")

View(cbind(names,coef(reslg1),true))
View(cbind(names,coef(reslg2),harder))

#really bad at sigma_uv no matter what
#if you have a starting value far from sigma_v it gets weird
#but sigma_v should actually be pretty easy to estimate


# manual check of values
test = pnorm(mtest)*(( 1/( 2*pi*(var(u)-cov(u,v)^2/var(v) ) ) )*
                       exp(-( 1/(2*(var(u) - cov(u,v)^2/var(v))) ) *
                             (log(y1) - y1*beta1 + x1*beta2 +alpha1*(y2 - x1*pi1 - x2*pi2))^2))
nI = I==FALSE
l.test = log(test)
l.test = (1-I)*l.test
l.test[is.infinite(l.test)] <- 0
l.test[is.nan(l.test)] <- 0

