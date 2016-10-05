##################################
# Replicating China 1 child paper
##################################

require(foreign)
require(plyr)

dat = read.dta("data_release2.dta")
attach(dat)

dat$post = ifelse(yr_of_b>1979, 1, 0)
dat$single = ifelse(only_child,1,0)
dat = dat[which(!is.na(yr_of_b)),]
yrb2 = ifelse(affected_78 == 1, NA, yr_of_b)
detach(dat); attach(dat)

# kicks out too many:
nona = dat[complete.cases(dat[ ,c(9:13)]),]

# replicate fig 1
# i need to kick out 5 ppl in 1975, 1 person in 1978 (not in 1978a), 1 person in 1983 
# % first born way, way off
ddply(dat,~yrb2,summarise,only=mean(only_child, na.rm = T)
      ,sib1=mean(siblings==1, na.rm = T)
      ,sib2=mean(siblings==2, na.rm = T)
      ,sib3=mean(siblings==3, na.rm = T)
      ,sib4=mean(siblings==4, na.rm = T)
      ,sibs=mean(siblings, na.rm = T)
      ,first=mean(first_born, na.rm = T)
      ,cuz=mean(cousins, na.rm = T)
      ,aunt=mean(aunts_uncles, na.rm = T)
      ,age=mean(age_days, na.rm = T)/365
      ,n = length(siblings)
      )

# need to change these to tobit (probit for the last one very different)
endog1 <- lm(dict_p_sent~post+male+uniplus+coll+bj+mcoll+muniplus)
endog2 <- lm(trust_p_sent~post+male+uniplus+coll+bj+mcoll+muniplus)
endog3 <- lm(perc_ret~post+male+uniplus+coll+bj+mcoll+muniplus)
endog4 <- lm(risk_invest~post+male+uniplus+coll+bj+mcoll+muniplus)
endog5 <- glm(compete~post+male+uniplus+coll+bj+mcoll+muniplus, family = binomial("probit"))


require(censReg)
reduced.form <- lm(single~post+male+uniplus+coll+bj+mcoll+muniplus)
summary(reduced.form)

consistent.tobit <- censReg(compete~fitted(reduced.form)+residuals(reduced.form))
summary(consistent.tobit)