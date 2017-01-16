source("wrapperFunc.R")
source("hurdleIVsim.R")
source("startFn.R")

dat = hurdle.IV.sim()

# input bad formula
hurdle.IV(formula = blob
          ,inst = inst1
          ,endog = endog
          ,exog = exog1
          ,data = dat
          ,endog_reg = list()
          ,start_val = F)

# more exogenous vars
dat = hurdle.IV.sim(formula = F,
                    pi = c(1,-1,.4,5),
                    gamma = c(-.2,.8,-.5,.07),
                    beta = c(.05,.06,.4,3),
                    exog_mean = c(1,1),
                    exog_sd = c(4,.6))

check = hurdle.IV(formula = y ~ endog + exog1 + exog2
          ,inst = inst1
          ,endog = endog
          ,exog = c(exog1,exog2)
          ,data = dat
          ,endog_reg = list()
          ,start_val = F)

# if you leave out endogenous regressor
check = hurdle.IV(formula = y ~  exog1 + exog2
                  ,inst = inst1
                  ,endog = endog
                  ,exog = c(exog1,exog2)
                  ,data = dat
                  ,endog_reg = list()
                  ,start_val = F)

# more insts
dat = hurdle.IV.sim(formula = F,
                    pi = c(1,-1,.4,5),
                    gamma = c(-.2,.8,.07),
                    beta = c(.05,.06,3),
                    z_mean = c(1,1),
                    z_sd = c(4,.6),
                    silent = F)

hurdle.IV(formula = y ~ endog + exog1 
          ,inst = c(inst1,inst2)
          ,endog = endog
          ,exog = exog1
          ,data = dat
          ,endog_reg = list()
          ,start_val = F)

# specify endogreg
dat = hurdle.IV.sim(formula = F,
                    pi = c(1,-1,.4,5),
                    gamma = c(-.2,.8,.07),
                    beta = c(.05,.06,3),
                    z_mean = c(1,1),
                    z_sd = c(4,.6),
                    silent = F)

hurdle.IV(formula = y ~ endog + exog1 
          ,inst = c(inst1,inst2)
          ,endog = endog
          ,exog = exog1
          ,data = dat
          ,endog_reg = list(endog~exog1 + inst1+inst2)
          ,start_val = F)

