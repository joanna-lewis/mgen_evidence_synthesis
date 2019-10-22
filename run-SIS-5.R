rm(list=ls())

set.seed(12345)
#setwd('/Users/Joanna/OneDrive - Imperial College London/backup/M_Genitalium/evidence_synthesis/stan')
#setwd("C:/Users/ucbpjle/OneDrive - Imperial College London/backup/M_Genitalium/evidence_synthesis/stan")

library(rstan)
options(mc.cores = parallel::detectCores()).
rstan_options(auto_write = TRUE)

dt <- list(
    prev_num_ct = 137, # Oakeshott 2010
    prev_denom_ct = 2377, # Oakeshott 2010
    
    prev_num_mg = 78, # Oakeshott 2010
    prev_denom_mg = 2378, # Oakeshott 2010
    
    N_persist = 1, # number of data points for persistence time
    persist_t = 16/12, # times for persistence follow-up
    persist_num = array(7), # numerators for persistence follow-up
    persist_denom = array(27), # denominators for persistence follow-up
    
    pid_given_infected_num_ct = 7, # Oakeshott 2010
		pid_given_infected_denom_ct = 70, # Oakeshott 2010
		pid_given_susceptible_num_ct = 31, # Oakeshott 2010
		pid_given_susceptible_denom_ct = 2114, # Oakeshott 2010 
		
		pid_given_infected_num_mg = 3, # Oakeshott 2010
		pid_given_infected_denom_mg = 77, # Oakeshott 2010
		pid_given_susceptible_num_mg = 36, # Oakeshott 2010
		pid_given_susceptible_denom_mg = 2169 # Oakeshott 2010
		
)

inits <- list(chain1 = list(NA), 
              chain2 = list(NA), 
              chain3 = list(NA), 
              chain4 = list(NA), 
              chain5 = list(NA)
              )

for(i in 1:5){
  inits[[i]] <- list(alpha_SC = exp(rnorm(1,0,1)), 
                     alpha_SM = exp(rnorm(1,0,1)), 
                     alpha_CS = exp(rnorm(1,0,1)), 
                     alpha_MS = exp(rnorm(1,0,1)), 
                     alpha_SP = exp(rnorm(1,0,1)), 
                     alpha_CP = exp(rnorm(1,0,1)), 
                     alpha_MP = exp(rnorm(1,0,1))
                     )
}

fitmgct <- stan('SIS-5.stan', 
              data = dt,
              chains = 5,
              iter = 20000,
              warmup = 5000,
              init = inits, #list(chain1=list(alpha_SI = 0.03, alpha_IS = 1, alpha_SP = 0.01, pid_diff = 0.05)),
              seed = 67890,
              control = list(max_treedepth = 10)
              )
	
opmgct <- extract(fitmgct)

plot(opmgct$lp__, type='l')
#pairs(fitmgct, pars=c("alpha_SC", "alpha_SM", "alpha_CS", "alpha_MS", "alpha_SP", "alpha_CP", "alpha_MP"))

par(mfrow=c(3,3), mar=c(5,4,1,1))

hist(opmgct$alpha_SC, n=20, freq=FALSE, main="", xlab=expression(alpha[SC]), xlim=c(0,0.1), ylim=c(0,80))
lines(seq(0,1,0.001), dgamma(seq(0,1,0.001),1,2), col="red")
quantile(opmgct$alpha_SC, p=c(0.5, 0.025, 0.975))

hist(opmgct$alpha_SM, n=20, freq=FALSE, main="", xlab=expression(alpha[SM]), xlim=c(0,0.1), ylim=c(0,80))
lines(seq(0,1,0.001), dgamma(seq(0,1,0.001),1,2), col="red")
quantile(opmgct$alpha_SM, p=c(0.5, 0.025, 0.975))

hc <- hist(opmgct$alpha_SC, n=50, plot=FALSE)
hm <- hist(opmgct$alpha_SM, n=50, plot=FALSE)

plot(rep(hc$breaks, each=2), c(0, rep(hc$density, each=2), 0), 
     type='l', col=2, xlim=c(0,0.1), xlab="Force of infection", ylab = "Density")
lines(rep(hm$breaks, each=2), c(0, rep(hm$density, each=2), 0), type='l', col=3)


hist(opmgct$alpha_CS, n=20, freq=FALSE, main="", xlab=expression(alpha[CS]), xlim=c(0,2.5), ylim=c(0,5))
lines(seq(0,1,0.001), dnorm(seq(0,1,0.001), 0.74, (0.89-0.61)/3.919928), col="red")
quantile(opmgct$alpha_CS, p=c(0.5, 0.025, 0.975))

hist(opmgct$alpha_MS, n=20, freq=FALSE, main="", xlab=expression(alpha[MS]), xlim=c(0,2.5), ylim=c(0,5))
lines(seq(0,10,0.001), dgamma(seq(0,10,0.001),1,2), col="red")
quantile(opmgct$alpha_MS, p=c(0.5, 0.025, 0.975))

hc <- hist(opmgct$alpha_CS, n=50, plot=FALSE)
hm <- hist(opmgct$alpha_MS, n=50, plot=FALSE)

plot(rep(hc$breaks, each=2), c(0, rep(hc$density, each=2), 0), 
     type='l', col=2, xlim=c(0,2), xlab="Clearance rate", ylab = "Density")
lines(rep(hm$breaks, each=2), c(0, rep(hm$density, each=2), 0), type='l', col=3)


hist(opmgct$alpha_SP, n=20, freq=FALSE, main="", xlab=expression(alpha[SP]), xlim=c(0,0.4), ylim=c(0,10))
lines(seq(0,1,0.001), dgamma(seq(0,1,0.001),1,2), col="red")
quantile(opmgct$alpha_SP, p=c(0.5, 0.025, 0.975))

hist(opmgct$alpha_CP, n=20, freq=FALSE, main = "", xlab=expression(alpha[CP]), xlim=c(0,0.4), ylim=c(0,10))
lines(seq(0,1,0.001), dgamma(seq(0,1,0.001),1,2), col="red")
quantile(opmgct$alpha_CP, p=c(0.5, 0.025, 0.975))

hist(opmgct$alpha_MP, n=20, freq=FALSE, main="", xlab=expression(alpha[MP]), xlim=c(0,0.4), ylim=c(0,10))
lines(seq(0,1,0.001), dgamma(seq(0,1,0.001),1,2), col="red")
quantile(opmgct$alpha_MP, p=c(0.5, 0.025, 0.975))

hs <- hist(opmgct$alpha_SP, n=20, plot=FALSE)
hc <- hist(opmgct$alpha_CP, n=100, plot=FALSE)
hm <- hist(opmgct$alpha_MP, n=100, plot=FALSE)

plot(rep(hs$breaks, each=2), c(0, rep(hs$density, each=2), 0), type='l', xlim=c(0,0.4), ylim=c(0,12))
lines(rep(hc$breaks, each=2), c(0, rep(hc$density, each=2), 0), type='l', col=2)
lines(rep(hm$breaks, each=2), c(0, rep(hm$density, each=2), 0), type='l', col=3)

quantile(opmgct$alpha_SP, p=c(0.5, 0.025,0.975))
quantile(opmgct$alpha_CP, p=c(0.5, 0.025,0.975))
quantile(opmgct$alpha_MP, p=c(0.5, 0.025,0.975))

# #quartz(width = 9, height = 6)
# par(mfrow=c(2,3))
# 
# hist(op0$alpha_SI, freq=FALSE)
# lines(seq(0,1,0.01), dlnorm(seq(0,1,0.01), 100, 10), col=2)
# 
# hist(op0$alpha_IS, freq=FALSE)
# lines(seq(0,10,0.01), dlnorm(seq(0,10,0.01), 100, 10), col=2)
# 
# hist(op$alpha_SP, freq=FALSE)
# lines(seq(0,1,0.01), dlnorm(seq(0,1,0.01), 100, 10), col=2)
# 
# hist(op0$pid_diff, freq=FALSE)
# lines(seq(0,1,0.01), dlnorm(seq(0,1,0.01), 100, 10), col=2)
# 
# hist(op$alpha_IP_mg, freq=FALSE)
# 
# # hist(op0$alpha_PI, freq=FALSE)
# # lines(seq(0,10,0.01), dgamma(seq(0,10,0.01), 0.001, 0.1), col=2)
# 
# hist(op0$prop_mgen, freq=FALSE)
# lines(seq(0,1,0.01), dbeta(seq(0,1,0.01), 2, 10), col=2)
# 
# hist(op0$alpha_PS, freq=FALSE)
# lines(seq(0,100,0.01), dgamma(seq(0,100,0.01), 50, 1), col=2)
# 
# hist(op0$alpha_IP /( op0$alpha_IP + op0$alpha_IS), freq=FALSE, xlim=c(0,1))
# quantile(op0$alpha_IP /( op0$alpha_IP + op0$alpha_IS))
# 
# quantile(op0$prev_num_sim, p=c(0.5, 0.025, 0.975)) # data: 29
# quantile(op0$pid_cases_sim, p=c(0.5, 0.025, 0.975)) # data: 40000
# 
# ######################################
# # look at expected PID, given different Mg prevalence
# ######################################
# 
# dt_plt <- array(dim=c(nrow=length(seq(0, 0.025, 0.0001)), 
#                       ncol=length(opmgct$alpha_IP_mg),
#                       4 # PID due to Mg, Ct, other causes; prob(reported)
#                       )) 
# 
# prev_mg <- seq(0, 0.025, 0.0001)  
# for(i in 1:nrow(dt_plt)){
#   dt_plt[i,,1] <- 100000*prev_mg[i]*op$alpha_IP_mg
#   dt_plt[i,,2] <- 100000*0.015*op$alpha_IP_ct    #1.5 % for 16-44s
#   dt_plt[i,,3] <- 100000*op$alpha_SP
# }
# dt_plt[,,4] <- matrix(rep(
#   rbeta(length(op$alpha_IP_mg), 12+9+1, 15+23 - (12+9) + 1),
#   each = nrow(dt_plt)
# ), nrow=length(seq(0, 0.025, 0.0001)), ncol=length(op$alpha_IP_mg))
#   
# plot(100*prev_mg, apply(dt_plt[,,1]*dt_plt[,,4], 1, median), type='l', ylim=c(0,2000),
#      xlab='Mgen prevalence (%)', ylab='PID per 100000 women per year')
# lines(100*prev_mg, apply(dt_plt[,,1]*dt_plt[,,4], 1, quantile, p=0.025), lty=2)
# lines(100*prev_mg, apply(dt_plt[,,1]*dt_plt[,,4], 1, quantile, p=0.975), lty=2)
# 
# lines(100*prev_mg, apply(apply(dt_plt[,,1:2],c(1,2),sum)*dt_plt[,,4], 1, quantile, p=0.5), col=2)
# lines(100*prev_mg, apply(apply(dt_plt[,,1:2],c(1,2),sum)*dt_plt[,,4], 1, quantile, p=0.025), lty=2, col=2)
# lines(100*prev_mg, apply(apply(dt_plt[,,1:2],c(1,2),sum)*dt_plt[,,4], 1, quantile, p=0.975), lty=2, col=2)
# 
# lines(100*prev_mg, apply(apply(dt_plt[,,1:3],c(1,2),sum)*dt_plt[,,4], 1, quantile, p=0.5), col=3)
# lines(100*prev_mg, apply(apply(dt_plt[,,1:3],c(1,2),sum)*dt_plt[,,4], 1, quantile, p=0.025), lty=2, col=3)
# lines(100*prev_mg, apply(apply(dt_plt[,,1:3],c(1,2),sum)*dt_plt[,,4], 1, quantile, p=0.975), lty=2, col=3)
# 
# 
# # add Natsal prevalence in 16-44-year-old women
# abline(v=1.3)
# polygon(c(0.9,0.9,1.9,1.9), c(-100, 2000, 2000, -100), border=NA, col=rgb(0,0,0,0.2))
# 
# #abline(v=1.7, lwd=2) # 16-24-year-old women
# 
# abline(h=304)
# polygon(c(-10,-10,10,10), c(289, 321, 321, 289), border=NA, col=rgb(0,0,0,0.2))
# #abline(h=289)
# #abline(h=321)
# 
# ###################
# # parameter posteriors
# ###################
# quantile(op0$alpha_SI, p=c(0.5, 0.025, 0.975)) 
# quantile(op0$alpha_SP, p=c(0.5, 0.025, 0.975)) 
# quantile(op0$alpha_IS, p=c(0.5, 0.025, 0.975))
# quantile(op0$alpha_IP, p=c(0.5, 0.025, 0.975))
# 
# ###################
# # prevalence
# ###################
# quantile(op0$I_star, p=c(0.5, 0.025, 0.975))
# 
###################
# proportion of Mg, Ct infections which result in PID
###################

quantile(opmgct$alpha_MP/(opmgct$alpha_MP + opmgct$alpha_MS), p=c(0.5, 0.025, 0.975))
quantile(opmgct$alpha_CP/(opmgct$alpha_CP + opmgct$alpha_CS), p=c(0.5, 0.025, 0.975))

quantile((opmgct$alpha_MP*(opmgct$alpha_CP + opmgct$alpha_CS)/
            (opmgct$alpha_CP*(opmgct$alpha_MP + opmgct$alpha_MS))), p=c(0.5, 0.025, 0.975))

quantile(
  (opmgct$alpha_MP/(opmgct$alpha_MP + opmgct$alpha_MS)) -
    (opmgct$alpha_CP/(opmgct$alpha_CP + opmgct$alpha_CS)), 
  p=c(0.5, 0.025, 0.975))

hm <- hist(100*opmgct$alpha_MP/(opmgct$alpha_MP + opmgct$alpha_MS), breaks=seq(0,100,0.5), plot=FALSE)
hc <- hist(100*opmgct$alpha_CP/(opmgct$alpha_CP + opmgct$alpha_CS), breaks=seq(0,100,0.5), plot=FALSE)

par(mfrow=c(1,2))
plot(rep(hm$breaks, each=2), c(0, rep(hm$density, each=2), 0), type="l", lwd = 5,
     main = "(a)",
     xlim = c(0,35), xlab = "Percentage of infections leading to PID",
     ylim=c(0,0.12), ylab = "Density")
lines(rep(hc$breaks, each=2), c(0, rep(hc$density, each=2), 0), col=2, lwd=3)
legend('topright', lwd=c(5,3), col=c(1,2), legend = c("Mgen", "Ct"), bty="n")

###################
# proportion of PID caused by M gen
###################
quantile((opmgct$M_star * opmgct$alpha_MP)/
            (opmgct$M_star*opmgct$alpha_MP + opmgct$C_star*opmgct$alpha_CP + opmgct$alpha_SP), 
        p=c(0.5, 0.025, 0.975))

quantile((opmgct$C_star * opmgct$alpha_CP)/
           (opmgct$M_star*opmgct$alpha_MP + opmgct$C_star*opmgct$alpha_CP + opmgct$alpha_SP), 
         p=c(0.5, 0.025, 0.975))

quantile(opmgct$alpha_SP/
           (opmgct$M_star*opmgct$alpha_MP + opmgct$C_star*opmgct$alpha_CP + opmgct$alpha_SP), 
         p=c(0.5, 0.025, 0.975))

hm <- hist(100*(opmgct$M_star * opmgct$alpha_MP)/
             (opmgct$M_star*opmgct$alpha_MP + opmgct$C_star*opmgct$alpha_CP + opmgct$alpha_SP),
           n=50, plot=FALSE)
hc <- hist(100*(opmgct$C_star * opmgct$alpha_CP)/
             (opmgct$M_star*opmgct$alpha_MP + opmgct$C_star*opmgct$alpha_CP + opmgct$alpha_SP),
           n=50, plot=FALSE)
hs <- hist(100*opmgct$alpha_SP/
             (opmgct$M_star*opmgct$alpha_MP + opmgct$C_star*opmgct$alpha_CP + opmgct$alpha_SP),
           n=50, plot=FALSE)

plot(rep(hm$breaks, each=2), c(0, rep(hm$density, each=2), 0), type="l", lwd=5,
     main = "(b)",
     xlim = c(0,100), xlab = "Percentage of PID attributable",
     ylim=c(0,0.07), ylab = "Density")
lines(rep(hc$breaks, each=2), c(0, rep(hc$density, each=2), 0), col=2, lwd=3)
lines(rep(hs$breaks, each=2), c(0, rep(hs$density, each=2), 0), col=3, lwd=1)
legend('topright', lwd=c(5,3,1), col=c(1,2,3), legend = c("Mgen", "Ct", "Other"), bty="n")

 
###################
# annual M gen PID cases (per person)
###################
100000*quantile(opmgct$M_star * opmgct$alpha_MP, p=c(0.5, 0.025, 0.975))
100000*quantile(opmgct$C_star * opmgct$alpha_CP, p=c(0.5, 0.025, 0.975))

#quantile(3432604.624 * op$I_star_mg * op$alpha_MP, p=c(0.5, 0.025, 0.975))

###################
# PID rate
###################
pid_rate <- apply(opmgct$natsal_pid, 2, quantile, p=c(0.5, 0.025, 0.975))
pid_rate_obs <- apply(opmgct$natsal_pid_obs, 2, quantile, p=c(0.5, 0.025, 0.975))

pid_mg_rate <- apply(opmgct$natsal_pid_mg, 2, quantile, p=c(0.5, 0.025, 0.975))
pid_ct_rate <- apply(opmgct$natsal_pid_ct, 2, quantile, p=c(0.5, 0.025, 0.975))

dt_price_16 <- data.frame(
  age = c("16-19", "20-24", "25-34", "35-44"),
  HES = c(1233, 3101, 9756, 10526),
  GPRD = c(5083, 8842, 14932, 9609),
  GUM = c(3212, 4399, 3919, 1388),
  pop = c(1199600, 1519100, 3502100, 3795600)*c(0.679, 0.944, 0.991, 0.991) # adjust for Natsal-2 proportion sexually active
)
dt_price_16$total_max <- dt_price_16$HES + dt_price_16$GPRD + dt_price_16$GUM
dt_price_16$rate_max <- 100000*dt_price_16$total_max/dt_price_16$pop
dt_price_16$total_min <- apply(dt_price_16[,c("HES","GPRD")], 1, max) + dt_price_16$GUM
dt_price_16$rate_min <- 100000*dt_price_16$total_min/dt_price_16$pop

#plot(1:4, pid_rate[1,]* (12+9)/(15+23), ylim=c(0, 1500), xlab="Age group", ylab = "Observed rate per 100000 women")
#arrows(x0=1:4, y0=pid_rate[2,]* (12+9)/(15+23), y1=pid_rate[3,]* (12+9)/(15+23), code=3, angle=90, length=0.05)

quartz(width=6, height=6)
par(mfrow=c(1,1), mar = c(5,5,4,2))

plot(1:4, 
     pid_rate_obs[1,], 
     xlim=c(0.5, 4.5), ylim=c(0, 1500), 
     xlab="Age group", ylab = "Observed PID rate per 100000 sexually-active women",
     xaxt = "n"
     )
axis(1, at=1:4, labels=dt_price_16$age)
arrows(x0=1:4, y0=pid_rate_obs[2,], y1=pid_rate_obs[3,], code=3, angle=90, length=0.05)

h1 <- hist(opmgct$natsal_pid_obs[,1], n=100, plot=FALSE)
h2 <- hist(opmgct$natsal_pid_obs[,2], n=100, plot=FALSE)
h3 <- hist(opmgct$natsal_pid_obs[,3], n=100, plot=FALSE)
h4 <- hist(opmgct$natsal_pid_obs[,4], n=100, plot=FALSE)

polygon(1 + c(0, 100*rep(h1$density, each=2), 0),rep(h1$breaks, each=2), col=rgb(0,0,0,0.3), border=NA)
polygon(1 - c(0, 100*rep(h1$density, each=2), 0),rep(h1$breaks, each=2), col=rgb(0,0,0,0.3), border=NA)
polygon(2 + c(0, 100*rep(h2$density, each=2), 0),rep(h2$breaks, each=2), col=rgb(0,0,0,0.3), border=NA)
polygon(2 - c(0, 100*rep(h2$density, each=2), 0),rep(h2$breaks, each=2), col=rgb(0,0,0,0.3), border=NA)
polygon(3 + c(0, 100*rep(h3$density, each=2), 0),rep(h3$breaks, each=2), col=rgb(0,0,0,0.3), border=NA)
polygon(3 - c(0, 100*rep(h3$density, each=2), 0),rep(h3$breaks, each=2), col=rgb(0,0,0,0.3), border=NA)
polygon(4 + c(0, 100*rep(h4$density, each=2), 0),rep(h4$breaks, each=2), col=rgb(0,0,0,0.3), border=NA)
polygon(4 - c(0, 100*rep(h4$density, each=2), 0),rep(h4$breaks, each=2), col=rgb(0,0,0,0.3), border=NA)

lines(1:4, dt_price_16$rate_max, lty=2)
lines(1:4, dt_price_16$rate_min, lty=2)

pid_ct_obs <- apply(opmgct$natsal_pid_ct_obs, 2, quantile, p=c(0.5, 0.025, 0.975))
points((1:4)-0.1, pid_ct_obs[1,], pch=2)
arrows((1:4)-0.1, y0=pid_ct_obs[2,], y1=pid_ct_obs[3,], code=3, length=0.05, angle=90)

h1c <- hist(opmgct$natsal_pid_ct_obs[,1], n=100, plot=FALSE)
h2c <- hist(opmgct$natsal_pid_ct_obs[,2], n=100, plot=FALSE)
h3c <- hist(opmgct$natsal_pid_ct_obs[,3], n=100, plot=FALSE)
h4c <- hist(opmgct$natsal_pid_ct_obs[,4], n=100, plot=FALSE)

# polygon(0.9 + 0.5*c(0, 100*rep(h1c$density, each=2), 0),rep(h1c$breaks, each=2), col=rgb(0,0,0,0.3), border=NA)
# polygon(0.9 - 0.5*c(0, 100*rep(h1c$density, each=2), 0),rep(h1c$breaks, each=2), col=rgb(0,0,0,0.3), border=NA)
# polygon(1.9 + 0.5*c(0, 100*rep(h2c$density, each=2), 0),rep(h2c$breaks, each=2), col=rgb(0,0,0,0.3), border=NA)
# polygon(1.9 - 0.5*c(0, 100*rep(h2c$density, each=2), 0),rep(h2c$breaks, each=2), col=rgb(0,0,0,0.3), border=NA)
# polygon(2.9 + 0.5*c(0, 100*rep(h3c$density, each=2), 0),rep(h3c$breaks, each=2), col=rgb(0,0,0,0.3), border=NA)
# polygon(2.9 - 0.5*c(0, 100*rep(h3c$density, each=2), 0),rep(h3c$breaks, each=2), col=rgb(0,0,0,0.3), border=NA)
# polygon(3.9 + 0.5*c(0, 100*rep(h4c$density, each=2), 0),rep(h4c$breaks, each=2), col=rgb(0,0,0,0.3), border=NA)
# polygon(3.9 - 0.5*c(0, 100*rep(h4c$density, each=2), 0),rep(h4c$breaks, each=2), col=rgb(0,0,0,0.3), border=NA)

pid_mg_obs <- apply(opmgct$natsal_pid_mg_obs, 2, quantile, p=c(0.5, 0.025, 0.975))
points((1:4)+0.1, pid_mg_obs[1,], pch=6)
arrows((1:4)+0.1, y0=pid_mg_obs[2,], y1=pid_mg_obs[3,], code=3, length=0.05, angle=90)

h1m <- hist(opmgct$natsal_pid_mg_obs[,1], n=100, plot=FALSE)
h2m <- hist(opmgct$natsal_pid_mg_obs[,2], n=100, plot=FALSE)
h3m <- hist(opmgct$natsal_pid_mg_obs[,3], n=100, plot=FALSE)
h4m <- hist(opmgct$natsal_pid_mg_obs[,4], n=100, plot=FALSE)

# polygon(1.1 + 0.5*c(0, 100*rep(h1m$density, each=2), 0),rep(h1m$breaks, each=2), col=rgb(0,0,0,0.3), border=NA)
# polygon(1.1 - 0.5*c(0, 100*rep(h1m$density, each=2), 0),rep(h1m$breaks, each=2), col=rgb(0,0,0,0.3), border=NA)
# polygon(2.1 + 0.5*c(0, 100*rep(h2m$density, each=2), 0),rep(h2m$breaks, each=2), col=rgb(0,0,0,0.3), border=NA)
# polygon(2.1 - 0.5*c(0, 100*rep(h2m$density, each=2), 0),rep(h2m$breaks, each=2), col=rgb(0,0,0,0.3), border=NA)
# polygon(3.1 + 0.5*c(0, 100*rep(h3m$density, each=2), 0),rep(h3m$breaks, each=2), col=rgb(0,0,0,0.3), border=NA)
# polygon(3.1 - 0.5*c(0, 100*rep(h3m$density, each=2), 0),rep(h3m$breaks, each=2), col=rgb(0,0,0,0.3), border=NA)
# polygon(4.1 + 0.5*c(0, 100*rep(h4m$density, each=2), 0),rep(h4m$breaks, each=2), col=rgb(0,0,0,0.3), border=NA)
# polygon(4.1 - 0.5*c(0, 100*rep(h4m$density, each=2), 0),rep(h4m$breaks, each=2), col=rgb(0,0,0,0.3), border=NA)

legend('topright', pch=c(1,2,6,NA), lty=c(0,0,0,2), cex=0.7,
       legend = c('Predicted: all causes', 'Predicted: Ct', 'Predicted: Mgen', 'Data: all causes'))

mtext("Predicted: all causes", side=3, at=-0.55, adj=0, line=2.3, cex=0.7)
axis(side=3, at=1:4, labels = paste(
  round(pid_rate_obs[1,]), 
  ' (', 
  round(pid_rate_obs[2,]), 
  '-', 
  round(pid_rate_obs[3,]),
  ')',
  sep = ""),
  tick = FALSE,
  line=1.3, cex.axis=0.7
)

mtext("Predicted: Ct", side=3, at=-0.55, adj=0, line=1.5, cex=0.7)
axis(side=3, at=1:4, labels = paste(
  round(pid_ct_obs[1,]), 
  ' (', 
  round(pid_ct_obs[2,]), 
  '-', 
  round(pid_ct_obs[3,]),
  ')',
  sep = ""),
  tick = FALSE,
  line=0.5, cex.axis=0.7
)

mtext("Predicted: Mgen", side=3, at=-0.55, adj=0, line=0.7, cex=0.7)
axis(side=3, at=1:4, labels = paste(
  round(pid_mg_obs[1,]), 
  ' (', 
  round(pid_mg_obs[2,]), 
  '-', 
  round(pid_mg_obs[3,]),
  ')',
  sep = ""),
  tick = FALSE,
  line=-0.3, cex.axis=0.7
)


# multi-panel plot with all posteriors
par(mfrow=c(3,3), mar=c(5,4,1,1))

hist(opmgct$alpha_SC, breaks = seq(0, 1, 0.002),
     freq=FALSE, main="", xlab=expression(alpha[SC]), xlim=c(0,0.1), ylim=c(0,80))
lines(seq(0,1,0.001), dgamma(seq(0,1,0.001),1,2), col="red")

hist(opmgct$alpha_SM, breaks = seq(0, 1, 0.002), 
     freq=FALSE, main="", xlab=expression(alpha[SM]), xlim=c(0,0.1), ylim=c(0,80))
lines(seq(0,1,0.001), dgamma(seq(0,1,0.001),1,2), col="red")

plot.new()

hist(opmgct$alpha_CS, breaks = seq(0, 3, 0.05), freq=FALSE, 
     main="", xlab=expression(alpha[CS]), xlim=c(0,2.5), ylim=c(0,5))
lines(seq(0,3,0.001), dnorm(seq(0,3,0.001), 0.74, (0.89-0.61)/3.919928), col="red")

hist(opmgct$alpha_MS, breaks = seq(0, 4, 0.05), 
     freq=FALSE, main="", xlab=expression(alpha[MS]), xlim=c(0,2.5), ylim=c(0,5))
lines(seq(0,10,0.001), dgamma(seq(0,10,0.001),1,2), col="red")

plot.new()

hist(opmgct$alpha_SP, breaks = seq(0, 1, 0.008), 
     freq=FALSE, main="", xlab=expression(alpha[SP]), xlim=c(0,0.4), ylim=c(0,10))
lines(seq(0,1,0.001), dgamma(seq(0,1,0.001),1,2), col="red")

hist(opmgct$alpha_CP, breaks = seq(0, 1, 0.008), 
     freq=FALSE, main = "", xlab=expression(alpha[CP]), xlim=c(0,0.4), ylim=c(0,10))
lines(seq(0,1,0.001), dgamma(seq(0,1,0.001),1,2), col="red")

hist(opmgct$alpha_MP, breaks = seq(0, 1, 0.008),
     freq=FALSE, main="", xlab=expression(alpha[MP]), xlim=c(0,0.4), ylim=c(0,10))
lines(seq(0,1,0.001), dgamma(seq(0,1,0.001),1,2), col="red")
