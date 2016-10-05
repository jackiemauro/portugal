############################################################################
#
# Store parameter values here to keep things consistent
#
############################################################################

lgiv_params = list(mux1 = c(1.3), sigx1 = c(sqrt(.2)),
              muz = 1.5,sigz = sqrt(.5),
              beta1 = c(-.05,-.15), gamma1 = c(.3,-.2), pi1 = c(.1,.2),
              beta2 = -.5, gamma2 = -.05, pi2 = .7,
              sig_v = .8, sig_u=.75, tau0=.04, tau1 = .02)

cr_params = list(mux1 = c(.2), sigx1 = c(sqrt(.2)),
                 muz = 1.5,sigz = sqrt(.5),
                 beta1 = c(-.05,-.5), gamma1 = c(-.5,-.2), pi1 = c(2),
                 beta2 = -.8, gamma2 = -.5, pi2 = .7,
                  sig_v = 1.8, sig_u=1.3, tau0=.04, tau1 = .02)

cr_params_mult = list(mux1 = c(.2,-.04), sigx1 = c(sqrt(.2),sqrt(.02)),
                 muz = c(.5,.03,-.01),sigz = c(sqrt(.5),sqrt(.05),sqrt(.05)),
                 beta1 = c(-.05,-.5,-.02), gamma1 = c(-.03,-.2,-.01), 
                 pi1 = list(c(.1,2,1),c(.01,.2),c(.7,.8)),
                 beta2 = c(-.8, -.1,.03), gamma2 = c(-.5,-.04,.05), 
                 pi2 = list(c(.7,.03,-.02),c(-.2,.09,-.09),c(1,-.02,.03)),
                 sig_v = c(1.8,1.2,.12), sig_u=1.3, 
                 tau0=c(.04,.04,.04), tau1 = c(.02,.02,.02))

lgiv_params_mult = list(mux1 = c(2,.3), sigx1 = c(sqrt(.2),sqrt(.02)),
                      muz = c(.5,1.5,-.01),sigz = c(sqrt(.5),sqrt(.05),sqrt(.05)),
                      beta1 = c(-.05,-.5,-.02), gamma1 = c(.3,-.2,-.01),
                      pi1 = list(c(.1,2,1),c(.01,.2),c(.7,.8)),
                      beta2 = c(-.8, -.1,.03), gamma2 = c(-.05,-.04,.05), 
                      pi2 = list(c(.7,.03,-.02),c(-.2,.09,-.09),c(1,-.02,.03)),
                      sig_v = c(.5,.8,.12), sig_u=.75, 
                      tau0=c(.04,.06,.09), tau1 = c(.01,.02,.04))


tob_params = cr_params
tob_params_mult = cr_params_mult

# check look
# gen_cragg1(n = 10000, params = cr_params)
# mean(as.numeric(y1==0))
# plot(y1[1,]~x2[2,])
# cr_params$gamma2 = 3
# gen_cragg1(n = 10000, params = cr_params)
# mean(as.numeric(y1==0))
# plot(y1~x2)


#tobit assumption values
assmp_tau1 <- 0
assmp_sigu <- 1
sig_y0_x2 <- 1 - cr_params$tau0^2/cr_params$sig_v^2
assmp_b2 <- (cr_params$gamma2 - cr_params$tau0/cr_params$sig_v^2)/sig_y0_x2
assmp_b1 <- cr_params$gamma1/sig_y0_x2
assmp_g1 <- cr_params$beta1*sig_y0_x2
assmp_g2 <- cr_params$beta2*sig_y0_x2 + cr_params$tau0/cr_params$sig_v^2

assmp_params = cr_params
assmp_params$tau1 = assmp_tau1
assmp_params$sig_u = assmp_sigu
assmp_params$gamma1 = assmp_g1
assmp_params$gamma2 = assmp_g2
