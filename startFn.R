start.val <- function(formula
                      ,endog_reg
                      ,data){
  lin.form = update(formula, log(.)~.)
  Iy = data$outcome>0
  data$y0 = as.numeric(Iy)
  linear = lm(lin.form, data = data[Iy,])
  prob.form = update(formula, y0 ~ . ) 
  probit = glm(prob.form, family = binomial("probit"), data = data)
  
  betas = coef(linear)
  y_sd = sd(linear$residuals)
  gammas = coef(probit)
  rho = cov(linear$residuals, probit$residuals[Iy])
  
  endogs = lapply(endog_reg, function(x) lm(x, data = data))
  pis = lapply(endogs, function(x) coef(x))
  endog_sd = lapply(endogs, function(x) sd(x$residuals))
  tau0 = lapply(endogs, function(x) cov(probit$residuals, x$residuals))
  tau1 = lapply(endogs, function(x) cov(linear$residuals, x$residuals[Iy]))
  
  return(list(beta = betas,
              gamma = gammas,
              pi = pis,
              endog_sd = unlist(endog_sd),
              y_sd = y_sd,
              tau0 = tau0,
              tau1 = tau1,
              rho = rho))
}