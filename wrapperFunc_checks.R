source("wrapperFunc.R")
source("hurdleIVsim.R")
source("startFn.R")
source("loglikFunc.R")
source("loglikFunc_cragg.R")

dat = hurdle.IV.sim(formula = F,
                    pi = c(0,0,0),
                    gamma = c(0,0,0),
                    beta = c(0,0,0),
                    endog_reg = list(),
                    exog_mean = 1,
                    exog_sd = 0,
                    z_mean = 3,
                    z_sd = 1,
                    endog_sd = 1,
                    y_sd = 1,
                    rho = 0,
                    tau0 = 0,
                    tau1 = 0,
                    n = 1000,
                    silent = F)

hurdle.IV(formula = y~endog + exog1
          ,inst = inst1
          ,endog = endog
          ,exog = exog1
          ,data = dat
          ,endog_reg = list()
          ,start_val = list())

dat = hurdle.IV.sim(formula = F,
                    pi = c(1,-1,-1),
                    gamma = c(1,.5,-2),
                    beta = c(.05,.08,.05),
                    endog_reg = list(),
                    exog_mean = 1,
                    exog_sd = 1,
                    z_mean = 3,
                    z_sd = 1,
                    endog_sd = 5,
                    y_sd = 2,
                    rho = .2,
                    tau0 = .3,
                    tau1 = .1,
                    n = 1000,
                    silent = F)

# give it the right answers
out = hurdle.IV(formula = y~endog + exog1
          ,inst = inst1
          ,endog = endog
          ,exog = exog1
          ,data = dat
          ,endog_reg = list()
          ,start_val = list(pi = list(c(1,-1,-1)),
                            gamma = c(1,.5,-2),
                            beta = c(.05,.08,.05),
                            endog_reg = list(),
                            exog_mean = 1,
                            exog_sd = 1,
                            z_mean = 3,
                            z_sd = 1,
                            endog_sd = 5,
                            y_sd = 2,
                            rho = .2,
                            tau0 = .3,
                            tau1 = .1)
          )


# input bad formula
hurdle.IV(formula = blob
          ,inst = inst1
          ,endog = endog
          ,exog = exog1
          ,data = dat
          ,endog_reg = list()
          ,start_val = F)

# check it gets the right betas and gammas
dat = hurdle.IV.sim(pi =c(1,1,0)
                    ,gamma = c(1,1,0)
                    ,beta = c(1,1,0))
hurdle.IV(formula = y~endog + exog1
          ,inst = inst1
          ,endog = endog
          ,exog = exog1
          ,data = dat
          ,endog_reg = list()
          ,start_val = list())

# simplest
dat = hurdle.IV.sim()
hurdle.IV(formula = y~ exog1 +endog
          ,inst = inst1
          ,endog = endog
          ,exog = exog1
          ,data = dat
          ,endog_reg = list()
          ,start_val = list())

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
          ,start_val = list())

# opposite betas and gammas
dat = hurdle.IV.sim(formula = F,
                    pi = c(1,-1,.4,5),
                    gamma = c(1,1,1,-.07),
                    beta = c(1,-1,-1,3),
                    exog_mean = c(1,1),
                    exog_sd = c(4,.6))

check = hurdle.IV(formula = y ~ endog + exog1 + exog2
                  ,inst = inst1
                  ,endog = endog
                  ,exog = c(exog1,exog2)
                  ,data = dat
                  ,endog_reg = list()
                  ,start_val = list())

# if you leave out endogenous regressor
check = hurdle.IV(formula = y ~  exog1 + exog2
                  ,inst = inst1
                  ,endog = endog
                  ,exog = c(exog1,exog2)
                  ,data = dat
                  ,endog_reg = list()
                  ,start_val = F)

# more insts -- shouldn't work
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

#add endog reg
dummy = rnorm(length(dat$y))
dat = data.frame(dat,dummy)

hurdle.IV(formula = y ~ endog + exog1 + dummy
          ,inst = c(inst1,inst2)
          ,endog = c(endog,exog1)
          ,exog = dummy
          ,data = dat
          ,endog_reg = list(endog~exog1 + inst1+inst2, exog1~inst2)
          ,start_val = F)

# cragg
dat = hurdle.IV.sim(type = "cragg")
hurdle.IV(formula = y~ exog1 +endog
          ,inst = inst1
          ,endog = endog
          ,exog = exog1
          ,data = dat
          ,endog_reg = list()
          ,start_val = list()
          ,type = "cragg")
