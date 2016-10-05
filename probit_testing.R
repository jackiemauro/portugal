n = 1000
out = rep(NA,n)
for(ii in 1:n){
  SigE = matrix(c(1,0,.1,0,3,.2,.1,.2,.5), byrow = T, ncol = 3)
  A = matrix(c(1,0,.1,0,1,-.2,0,0,1), byrow = T, ncol = 3)
  Sig = t(A)%*%SigE%*%A
  
  euv = mvrnorm(n=100, mu = c(rep(0,dim(Sig)[1])), Sigma = Sig)
  eta = euv[,1]
  u = euv[,2]
  v = euv[,3]
  
  t01 = rnorm(100, .5, .1); t02 = rnorm(100, 1.5, .2)
  t1 = rnorm(100, 1.3, 3)
  t2 = .1 + 2*t1 + .7*t01 + .03*t02 + rnorm(100,0,.5)
  v2 = 2^2*var(t1) + .7^2*var(t01) + .03^2*var(t02) + .5^2
  
  t3 = -.03 + -.2*t1- .05*t2 
  v3 = (-.2)^2*var(t1) + (-.05)^2*var(t2)
  out[ii] = (v3 - var(t3))/var(t3)
}
hist(out)

t01 = rnorm(100, .5, .1); t02 = rnorm(100, 1.5, .2)
t1 = rnorm(100, 1.3, 3)
t2 = .1 + 2*t1 + .7*t01 + .03*t02 + rnorm(100,0,.5)
v2 = 2^2*var(t1) + .7^2*var(t01) + .03^2*var(t02) + .5^2
t3 = -.03 + -.2*t1- .05*t2 
v3 = (-.2)^2*var(t1) + (-.04)^2*v2
pr = pnorm(t3)

y = rbinom(100,1,pr)
df = data.frame(y,t1,t2)

test = glm(y~t1 + t2, data = df, family = binomial(link = 'probit'))

############ simpler

t1 = rnorm(100, 1.3, 1)
t2 = rnorm(100, 3, 3)
vars = matrix(c(t1,t2), ncol = 2, byrow = F)
coefs = matrix(c(.2, -.5), ncol = 1)
t3 = 1 + vars%*%coefs + rnorm(100)
sig2 = c(var(t1),var(t2))%*%coefs^2 + 1
(var(t3) - sig2)/sig2
pr = pnorm(t3)

y = rbinom(100,1,pr)
df = data.frame(y,t1,t2)

test = glm(y~t1 + t2, data = df, family = binomial(link = 'probit'))

############ simplest
t1 = rnorm(100, 2, .1)
t2 = rnorm(100, 3, .3)
t3 = 3 + 6*t1 - 4*t2 + rnorm(100,0,1)
pr = pnorm(t3)

y = rbinom(100,1,pr)
df = data.frame(y,t1,t2)

test = glm(y~t1 + t2, data = df, family = binomial(link = 'probit'))