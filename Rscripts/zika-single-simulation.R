# - - - - - - - - - - - - - - - - - - - - - - - 
# Code to run one simulation of ZIKV transmission model
# Author: Alasdair Henderson
# github.com/
# - - - - - - - - - - - - - - - - - - - - - - - 

zika_single_sim <- function(transmission_rate){
zika_data <- data.frame(date=as.Date(date.vals), Zika2016=y.vals)
zika_data$date <- as.Date(zika_data$date)

# Change parameter values -------------------------------------------------
par(mfrow=c(2,1), mar = c(3,4,1,3))

theta <- c(thetatab[1,],thetaAlltab[1,iiH,], theta_denv)
theta_init <- theta_initAlltab[1,iiH,]


introductions_mid <- theta[["intro_mid"]]
introductions_width <- theta[["intro_width"]]
introductions_base <- theta[["intro_base"]]
## total introductions
total_intro <- 4 * introductions_width * introductions_base # Integral of 4*base*exp(-(time-mid)/width)/(1+exp(-(time-mid)/width))^2 over -infty/infty is 4*base*width

theta[["beta_h"]] <- transmission_rate
# set initial conditions
init1 <- c(
  s_init=theta_init[["s_init"]],e_init=theta_init[["i1_init"]],i_init=theta_init[["i1_init"]],r_init=theta_init[["r_init"]],c_init=0,
  sd_init=theta_init[["sd_init"]],ed_init=theta_init[["ed_init"]],id_init=theta_init[["id_init"]],t1d_init=theta_init[["t1d_init"]],t2d_init=theta_init[["t2d_init"]],cd_init=0,
  sm_init=theta_init[["sm_init"]],em_init=theta_init[["em_init"]],im_init=theta_init[["im_init"]])

# time series, date series and case data
# adjust when DENV and ZIKV introductions begin (depending on when model start date is): n.b. DENV3 fixed to 2013-10-27
theta[["denv_start_point"]] <- as.Date("2013-10-27") - startdate
theta[["zika_start_point"]] <- theta[["intro_mid"]]
if(!is.na(theta[['epsilon']])){
  epsilon <- theta[['epsilon']]}else{
    epsilon <- 0}
if(!is.na(theta[['rho']])){
  theta[["rho"]] <- 1/theta[['rho']]}else{
    theta[["rho"]] <- 0}
if(!is.na(theta[['omega_d']])){
  theta[["omega_d"]] <- 1/theta[['omega_d']]}else{
    theta[["omega_d"]] <- 0}
if(!is.na(theta[['chi']])){
  theta[["chi"]] <- theta[['chi']]}else{
    theta[["chi"]] <- 0}
# Run simulation ----------------------------------------------------------
output <- simulate_deterministic_noage_DENVimm(theta, init1, time.vals.sim=time.vals)

# Match compartment states at sim.vals time
S_traj <- output[match(time.vals,output$time),"s_init"]
X_traj <- output[match(time.vals,output$time),"sm_init"]
R_traj <- output[match(time.vals,output$time),"r_init"]
I_traj <- output[match(time.vals,output$time),"i_init"]
cases1 <- output[match(time.vals,output$time),"c_init"]
casesD <- output[match(time.vals,output$time),"cd_init"]
SD <- output[match(time.vals,output$time),"sd_init"]
ED <- output[match(time.vals,output$time),"ed_init"]
ID <- output[match(time.vals,output$time),"id_init"]
RD <- output[match(time.vals,output$time),"rd_init"]
casecountD <- casesD-c(0,casesD[1:(length(time.vals)-1)])
casecount <- cases1-c(0,cases1[1:(length(time.vals)-1)])
casecount[casecount<0] <- 0
cases_est <- ReportC(casecount, theta[["rep"]], 1/theta[["repvol"]] )
cases_est[cases_est==0] <- NA

# calculate R0 ------------------------------------------------------------
tMax <- length(casecount)
t.start = 0
time.V = (1:tMax)*7  
model_st <- startdate
date0 = (model_st-date.vals[1]) %>% as.numeric() 

beta_ii <- seasonal_f(time.V[1:tMax], date0, amp=theta[["beta_v_amp"]], mid=theta[["beta_v_mid"]])
#decline_ii <- decline_f(time.V[1:tMax], mid=theta[["beta_mid"]], width=theta[["beta_grad"]], base=theta[["beta_base"]])
decline_ii <- control_f(time.V[1:tMax], base=theta[["beta_base"]], grad=theta[["beta_grad"]], mid=theta[["beta_mid"]], mid2=theta[["beta_mid"]]+theta[["beta_width"]], width=theta[["beta_width"]])
b_vary = beta_ii#*decline_ii

s_pick = S_traj[1:tMax]/342000
x_pick = X_traj[1:tMax] 
c_pick = casecount[1:tMax]/342000
r_pick = R_traj[1:tMax]/342000
theta[["Inf."]] <- theta[["Inf"]]
output_rr <- calculate_r0(th_in=as.data.frame(t(theta)),sus_c=s_pick,sus_a=0,sm_c=x_pick,sm_a=0,b_vary=b_vary,control=decline_ii)

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
  if(date < min(date.vals)){ 
    seroP[i] <- epsilon
    binom.lik[i] <- (dbinom(nLUM[lum.y==sero.y[i]], size=nPOP[lum.y==sero.y[i]], prob=seroP[i], log = T))
  }else{ 
    seroP[i] <-  (min(R_traj[date.vals<date+3.5 & date.vals>date-3.5])/theta[["npop"]]) + 
      (1 - min(R_traj[date.vals<date+3.5 & date.vals>date-3.5])/theta[["npop"]])*epsilon
    binom.lik[i] <- (dbinom(nLUM[lum.y==sero.y[i]], size=nPOP[lum.y==sero.y[i]], prob=seroP[i], log = T))
  }
  i <- i+1
}

ln.full <- length(y.vals)
likelihood <- sum(binom.lik) + sum(log(dnbinom(y.vals,
                                               mu=theta[["rep"]]*(casecount),
                                               size=1/theta[["repvol"]]))) 

# Plot simulation ---------------------------------------------------------
plot(date.vals, casecount, type='l', col=2, xlab="Year", ylab="Infections (cases)",ylim=c(0,2e4),
     main=paste0("Start: ", startdate+theta[["intro_mid"]],"   |   lik:", signif(likelihood,4),"  |  medR0: ", signif(median(r0_post),3), " | intro:", round(total_intro,0)))
lines(date.vals,casecountD)
par(new=T)
plot(date.vals, y.vals, type='h', col=4, yaxt='n', xaxt='n', xlab="", ylab="", ylim=c(0,5))
lines(date.vals,sapply(time.vals, function(xx){1-intro_f(xx, mid = theta[["zika_start_point"]], width = theta[["intro_width"]], base = theta[["intro_base"]])}),
      type='l',col=alpha(9,0.4), yaxt='n', xaxt='n', xlab="", ylab="")
axis(side=4)
mtext("Introductions", side=4, cex = 0.7, padj=3)

# Plot seasonality
plot(date.vals,rr_post,type='l',col=1, ylim=c(0,3), xlab="", ylab="")
lines(date.vals,decline_ii, type='l',col=2, yaxt='n', xaxt='n')
par(new=T)
plot(date.vals, (R_traj/342000)+theta[["epsilon"]], type='l', col=22, yaxt='n', xaxt='n', xlab="", ylab="", ylim=c(0,1))
points(seroposdates, (nLUM/nPOP), col=4, pch=1)
axis(side=4)

return(likelihood)
}
