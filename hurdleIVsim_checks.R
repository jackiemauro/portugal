# testing hurdleIVsim file
source("hurdleIVsim.R")

#basic
check = hurdle.IV.sim()

# doesn't take formulae
check = hurdle.IV.sim(formula = "dog")

# diff dimensions
check = hurdle.IV.sim(n = 10)
dim(check)

# more exogenous variables
check = hurdle.IV.sim(formula = F,
                      pi = c(1,-1,.4,5),
                      gamma = c(-.2,.8,-.5,.07),
                      beta = c(.05,.06,.4,3),
                      endog_reg = F,
                      exog_mean = c(1,1),
                      exog_sd = c(4,.6),
                      z_mean = 3,
                      z_sd = 1,
                      endog_sd = 3,
                      y_sd = 5,
                      n = 100)

summary(lm(endog~exog1 + exog2 + inst1, data= check))
summary(glm(y0~exog1 + exog2 + x2, family = binomial("probit"),data= check))

# messing with variances
check = hurdle.IV.sim(formula = F,
                      pi = c(1,-1.5,1,5),
                      gamma = c(-.2,.8,-.5,.7),
                      beta = c(.05,.06,.4,.07),
                      endog_reg = F,
                      exog_mean = c(1,1),
                      exog_sd = c(1,1),
                      z_mean = 3,
                      z_sd = .5,
                      endog_sd = 3,
                      y_sd = 2,
                      rho = 0, 
                      tau0 = .4,
                      tau1 = -10,
                      n = 100)

# messing with vector lengths
check = hurdle.IV.sim(formula = F,
                      pi = c(1,-1.5,1,5),
                      gamma = c(-.2,.8,-.5,.7),
                      beta = c(.05,.06,.4,.07),
                      endog_reg = F,
                      exog_mean = c(1,1),
                      exog_sd = c(1,1),
                      z_mean = 3,
                      z_sd = .5,
                      endog_sd = 3,
                      y_sd = 2,
                      rho = 0, 
                      tau0 = c(.4,5),
                      tau1 = .7,
                      n = 100)

