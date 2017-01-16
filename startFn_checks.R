#startFn checks
source("startFn.R")
source("hurdleIVsim.R")

# basic
dat2 = hurdle.IV.sim()
start.val(y~endog+exog1
          , endog_reg = list(endog~exog1+inst1)
          , data = dat2)

# messing up zeros
dat2$y0 = sample(c(0,1),length(dat$y0), replace = T, prob = c(.3,.7))
start.val(y~endog+exog1
          , endog_reg = list(endog~exog1+inst1)
          , data = dat2)

# more than one endog reg works
