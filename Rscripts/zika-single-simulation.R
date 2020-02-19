# - - - - - - - - - - - - - - - - - - - - - - - 
# MCMC analysis of Zika transmission in Fiji  
# Author: Alasdair Henderson
# github.com/
# - - - - - - - - - - - - - - - - - - - - - - - 

# Load packages and set up pathogen -----------------------------------------------
library(deSolve)
library(truncnorm)
library(magrittr)
library(ggplot2)
#rm(list=ls())

# set working directory ---------------------------------------------------
if(Sys.info()["user"]=="lsh1510922"){
  setwd("~/Documents/fiji-arbovirus/Zika_transmission_model")
}

# load functions ----------------------------------------------------------
source("preamble_conditions_zikafit.R")
source("R_functions/simulate_deterministic_noage_denvimmunity.R")
source("R_functions/calculateR0.R")
source("R_functions/seasonal_decline.R")
#
source("prior_distributions.R")
Virus="Zika"
source("load_serology_data.R")
zika_data <- read.csv("data_sets/Central_Fiji_2016Z_timeseries.csv")

# Fixed start date of the model  ------------------------------------------
model_start_date <- as.Date("2013-01-01")

# Change parameter values -------------------------------------------------
par(mfrow=c(2,1), mar = c(3,4,1,3))
cross_protection <- 0.3
control_measures <- 0.86528859

transmission_rate <- 0.33
 introductions_mid <- 300
 introductions_width <- 10
 introductions_base <- 1

## total introductions
total_intro <- 4 * introductions_width * introductions_base # Integral of 4*base*exp(-(time-mid)/width)/(1+exp(-(time-mid)/width))^2 over -infty/infty is 4*base*width
total_intro

# Need to fix  

# Set up data and initial conditions --------------------------------------
theta_star <- c("npop"=342000)

## beta_base/grad/mid=effect/length/midpoint of control measures (logistic function)
## chi=effect of cross protection [0,1] | omega_d=length of temporary cross-protection
## rho=lenght of time antibodies are detectable by MIA | epsilon=sensitivity of assay (prob false +)
## beta_d/gamma_d/alpha_d=transmission rate/infectious period/incubation period of DENV3 outbreak
thetaA_star <- c("beta_h"=NA, "beta_v"=NA, "beta_v_amp"=0.38106713, "beta_v_mid"=0.82821561,
                 "rep"=0.001, "repvol"=1, "Vex"=1/15, "Exp"=1/6.1, "MuV"=1/8.1, "Inf"=1/5, 
                 "beta_base"=NA, "beta_grad"=0.2, "beta_mid"=470, "beta_width"=20.03805649,
                 "tau"=1, "m"=0.56,
                 "chi"=NA, "omega_d"=1/30, "rho"=1/400, "epsilon"=0, #0.067,
                 "beta_d"=0.28866247, "gamma_d"=0.09821938, "alpha_d"=0.41303842)

thetaA_star[["beta_h"]] <- transmission_rate
thetaA_star[["beta_v"]] <- 1
thetaA_star[["chi"]] <- cross_protection
thetaA_star[["beta_base"]] <- control_measures
thetaA_star[["intro_base"]] <- introductions_base
thetaA_star[["intro_width"]] <- introductions_width
theta <- c(theta_star, thetaA_star)
# initial conditions
init1=c(
  s_init=theta[["npop"]], e_init=0, i_init=0, r_init=0, c_init=0,
  sd_init=228635.2028, ed_init=0, id_init=162.7972, t1d_init=0, t2d_init=0, cd_init=0,
  sm_init=1, em_init=0, im_init=0)

# time series, date series and case data
time_vals <- seq(7, as.numeric(as.Date("2018-01-08")-model_start_date), 7)
date_vals <- as.Date(time_vals, origin=model_start_date)
y_vals <- c(rep(0, (as.Date(zika_data$date[1]) - model_start_date)/7), zika_data$Zika2016)

# adjust when DENV and ZIKV introductions begin (depending on when model start date is): n.b. DENV3 fixed to 2013-10-27
theta[["denv_start_point"]] <- as.Date("2013-10-27") - model_start_date
theta[["zika_start_point"]] <- introductions_mid

#plot(date_vals, sapply(time_vals, function(xx){seasonal_f(xx, amp=theta[["beta_v_amp"]], mid=theta[["beta_v_mid"]])}), type="l")

# Run simulation ----------------------------------------------------------
output <- simulate_deterministic_noage_DENVimm(theta, init1, time.vals.sim=time_vals)


# Match compartment states at sim.vals time
S_traj <- output[match(time_vals,output$time),"s_init"]
X_traj <- output[match(time_vals,output$time),"sm_init"]
R_traj <- output[match(time_vals,output$time),"r_init"]
I_traj <- output[match(time_vals,output$time),"i_init"]
cases1 <- output[match(time_vals,output$time),"c_init"]
casesD <- output[match(time_vals,output$time),"cd_init"]
SD <- output[match(time_vals,output$time),"sd_init"]
ED <- output[match(time_vals,output$time),"ed_init"]
ID <- output[match(time_vals,output$time),"id_init"]
RD <- output[match(time_vals,output$time),"rd_init"]
casecountD <- casesD-c(0,casesD[1:(length(time_vals)-1)])
casecount <- cases1-c(0,cases1[1:(length(time_vals)-1)])
casecount[casecount<0] <- 0
cases_est <- ReportC(casecount, thetaA_star[["rep"]], 1/thetaA_star[["repvol"]] )
cases_est[cases_est==0] <- NA

# calculate R0 ------------------------------------------------------------
tMax <- length(casecount)
t.start = 0
time.V = (1:tMax)*7  
model_st <- model_start_date
date0 = (model_st-date_vals[1]) %>% as.numeric() 

beta_ii <- seasonal_f(time.V[1:tMax], date0, amp=thetaA_star[["beta_v_amp"]], mid=thetaA_star[["beta_v_mid"]])
#decline_ii <- decline_f(time.V[1:tMax], mid=thetaA_star[["beta_mid"]], width=thetaA_star[["beta_grad"]], base=thetaA_star[["beta_base"]])
decline_ii <- control_f(time.V[1:tMax], base=theta[["beta_base"]], grad=theta[["beta_grad"]], mid=theta[["beta_mid"]], mid2=theta[["beta_mid"]]+theta[["beta_width"]], width=theta[["beta_width"]])
b_vary = beta_ii#*decline_ii

s_pick = S_traj[1:tMax]/342000
x_pick = X_traj[1:tMax] 
c_pick = casecount[1:tMax]/342000
r_pick = R_traj[1:tMax]/342000
thetaA_star[["Inf."]] <- thetaA_star[["Inf"]]
output_rr = calculate_r0_adam(th_in=as.data.frame(t(thetaA_star)),sus_c=s_pick,sus_a=0,sm_c=x_pick,sm_a=0,b_vary=b_vary,control=decline_ii)

start.rr <- output_rr$rr_out[min(which(output_rr$rr_out>0))]
output_rr$rr_out[1:(min(which(output_rr$rr_out>0))-1)] <- start.rr
output_rr$r0_out[1:(min(which(output_rr$r0_out>0))-1)] <- start.rr
rr_post = output_rr$rr_out
rr_post[length(rr_post)] = rr_post[length(rr_post)-1]
r0_post = output_rr$r0_out
r0_post[length(r0_post)] = r0_post[length(r0_post)-1]
decline_post= decline_ii  

# Compute likelihood ------------------------------------------------------
# Calculate seropositivity at pre-specified dates and corresponding likelihood
epsilon <- theta[["epsilon"]]
i=1; seroP=NULL; binom.lik=NULL
sero.years <- format(as.Date(seroposdates, format="%d/%m/%Y"),"%Y")
sero.y <- substr(sero.years,3,4)
lum.y <- c("13","15","17")
for(date in seroposdates){
  if(date < min(date_vals)){ 
    seroP[i] <- epsilon
    binom.lik[i] <- (dbinom(nLUM[lum.y==sero.y[i]], size=nPOP[lum.y==sero.y[i]], prob=seroP[i], log = T))
  }else{ 
    seroP[i] <-  (min(R_traj[date_vals<date+3.5 & date_vals>date-3.5])/theta[["npop"]]) + 
      (1 - min(R_traj[date_vals<date+3.5 & date_vals>date-3.5])/theta[["npop"]])*epsilon
    binom.lik[i] <- (dbinom(nLUM[lum.y==sero.y[i]], size=nPOP[lum.y==sero.y[i]], prob=seroP[i], log = T))
  }
  i <- i+1
}

ln.full <- length(y_vals)
likelihood <- sum(binom.lik) + sum(log(dnbinom(y_vals,
                                               mu=theta[["rep"]]*(casecount),
                                               size=1/theta[["repvol"]]))) 

# likelihood <- sum(binom.lik) + sum(log(dpois(y_vals,
#                                     lambda=theta[["rep"]]*(casecount))))


# Plot simulation ---------------------------------------------------------
plot(date_vals, casecount, type='l', col=2, xlab="Year", ylab="Infections (cases)",ylim=c(0,2e4),
     main=paste0("Start: ", model_start_date+introductions_mid,"   |   lik:", signif(likelihood,4),"  |  medR0: ", signif(median(r0_post),3), " | intro:", round(total_intro,0)))
lines(date_vals,casecountD)
par(new=T)
plot(date_vals, y_vals, type='h', col=4, yaxt='n', xaxt='n', xlab="", ylab="", ylim=c(0,5))
#points(date_vals, cases_est, type='p', col=alpha(4,0.5))
#par(new=T)

#axis(side=4)
#abline(h=1, lty=2, col=alpha(9,0.1))
#par(new=T)
lines(date_vals,sapply(time_vals, function(xx){1-intro_f(xx, mid = theta[["zika_start_point"]], width = theta[["intro_width"]], base = theta[["intro_base"]])}),
      type='l',col=alpha(9,0.4), yaxt='n', xaxt='n', xlab="", ylab="")
axis(side=4)
mtext("Introductions", side=4, cex = 0.7, padj=3)
#print(likelihood)

# Plot seasonality
plot(date_vals,rr_post,type='l',col=alpha(9,0.8), ylim=c(0,3), xlab="", ylab="")
#lines(date_vals,r0_post, type='l',col=alpha(12,0.8), yaxt='n', xaxt='n')
lines(date_vals,decline_ii, type='l',col=alpha(2,0.8), yaxt='n', xaxt='n')
par(new=T)
plot(date_vals, (R_traj/342000)+theta[["epsilon"]], type='l', col=alpha(22,0.6), yaxt='n', xaxt='n', xlab="", ylab="", ylim=c(0,1))
points(seroposdates, (nLUM/nPOP), col=4, pch=1)
axis(side=4)

get_lik <- function(intro){
  theta[["zika_start_point"]] <- intro
  output <- simulate_deterministic_noage_DENVimm(theta, init1, time.vals.sim=time_vals)
  casecount <- cases1-c(0,cases1[1:(length(time_vals)-1)])
  casecount[casecount<0] <- 0
  R_traj <- output[match(time_vals,output$time),"r_init"]
  
  # Compute likelihood ------------------------------------------------------
  # Calculate seropositivity at pre-specified dates and corresponding likelihood
  epsilon <- theta[["epsilon"]]
  i=1; seroP=NULL; binom.lik=NULL
  sero.years <- format(as.Date(seroposdates, format="%d/%m/%Y"),"%Y")
  sero.y <- substr(sero.years,3,4)
  lum.y <- c("13","15","17")
  for(date in seroposdates){
    if(date < min(date_vals)){ 
      seroP[i] <- epsilon
      binom.lik[i] <- (dbinom(nLUM[lum.y==sero.y[i]], size=nPOP[lum.y==sero.y[i]], prob=seroP[i], log = T))
    }else{ 
      seroP[i] <-  (min(R_traj[date_vals<date+3.5 & date_vals>date-3.5])/theta[["npop"]]) + 
        (1 - min(R_traj[date_vals<date+3.5 & date_vals>date-3.5])/theta[["npop"]])*epsilon
      binom.lik[i] <- (dbinom(nLUM[lum.y==sero.y[i]], size=nPOP[lum.y==sero.y[i]], prob=seroP[i], log = T))
    }
    i <- i+1
  }
  
  ln.full <- length(y_vals)
  likelihood <- sum(binom.lik) + sum(log(dnbinom(y_vals,
                                                 mu=theta[["rep"]]*(casecount),
                                                 size=1/theta[["repvol"]])))
  likelihood
}

#likelihoods <- sapply(1:400, function(xx){get_lik(xx)})


#print(c(ii, likelihood))
total_intro
print(model_start_date+theta[["zika_start_point"]])
