# - - - - - - - - - - - - - - - - - - - - - - - 
# Code to run one simulation of ZIKV transmission model
# Author: Alasdair Henderson
# github.com/
# - - - - - - - - - - - - - - - - - - - - - - - 

denv3_single_sim <- function(transmission_rate){
  denv3_data <- data.frame(date=as.Date(date.vals), cases = y.vals)
  denv3_data$date <- as.Date(denv3_data$date)
  
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
    sd_init=theta_init[["sd_init"]],ed_init=theta_init[["ed_init"]],id_init=theta_init[["id_init"]],t1d_init=theta_init[["t1d_init"]],t2d_init=theta_init[["t2d_init"]],cd_init=0)
  
  # time series, date series and case data
  # adjust when DENV and ZIKV introductions begin (depending on when model start date is): n.b. DENV3 fixed to 2013-10-27
  theta[["denv_start_point"]] <- as.Date("2013-10-27") - startdate
  theta[["zika_start_point"]] <- theta[["intro_mid"]]
  if(!is.na(theta[['epsilon']])){
    epsilon <- theta[['epsilon']]}else{
      epsilon <- 0}
  if(!is.na(theta[['alpha']])){
    alpha <- theta[['alpha']]}else{
      alpha <- 1}
  if(!is.na(theta[['rho']])){
    theta[["rho"]] <- 1/theta[['rho']]}else{
      theta[["rho"]] <- 0}
  if(!is.na(theta[['omega_d']])){
    theta[["omega_d"]] <- 1/theta[['omega_d']]}else{
      theta[["omega_d"]] <- 0}
  if(!is.na(theta[['chi']])){
    theta[["chi"]] <- theta[['chi']]}else{
      theta[["chi"]] <- 0}
  if(!is.na(theta[['mu']])){
    theta[["mu"]] <- 1/(theta[['mu']]*365.25)
    theta[["eta"]] <- theta[['mu']]*2.5
  }else{
    theta[["eta"]] <- 0
    theta[["mu"]] <- 0}
  # Run simulation ----------------------------------------------------------
  output <- zikv_model_ode(theta, init1, time.vals.sim=time.vals)
  
  # Match compartment states at sim.vals time
  S_traj <- output[match(time.vals,output$time),"s_init"]
  R_traj <- output[match(time.vals,output$time),"r_init"]
  I_traj <- output[match(time.vals,output$time),"i_init"]
  cases1 <- output[match(time.vals,output$time),"c_init"]
  casesD <- output[match(time.vals,output$time),"cd_init"]
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
  decline_ii <- control_f(time.V[1:tMax], base=theta[["beta_base"]], mid=theta[["beta_mid"]], width=theta[["beta_width"]])
  b_vary = beta_ii
  
  s_pick = S_traj[1:tMax]/342000
  c_pick = casecount[1:tMax]/342000
  r_pick = R_traj[1:tMax]/342000
  theta[["Inf."]] <- theta[["Inf"]]
  theta[["mu"]] <- 67.25
  
  output_rr <- calculate_r0(th_in = as.data.frame(t(theta)),sus = s_pick, b_vary = b_vary)
  
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
  alpha <- theta[["alpha"]]

  i=1; seroP=NULL; binom.lik=NULL
  if(include.sero.likelihood==T){
    for(date in seroposdates){
      time_elapsed <- min(time.vals[date.vals<date+(dt/2) & date.vals>date-(dt/2)])
      adjusted_pop <- theta[["npop"]]-(time_elapsed*theta[["mu"]]*theta[["npop"]])+(time_elapsed*theta[["eta"]]*theta[["npop"]])
      modelled_R <- (min(R_traj[date.vals<date+(dt/2) & date.vals>date-(dt/2)])/adjusted_pop)
      detected_positives <- modelled_R*alpha    ## identify alpha% of the actual positives
      false_positives <- (1-modelled_R)*epsilon ## falsely record epsilon% of the actual negatives as positives
      seroP[i] <- detected_positives + false_positives
      seroP[i] <- (min(R_traj[date.vals<date+(dt/2) & date.vals>date-(dt/2)])/theta[["npop"]]) + 
        (1 - min(R_traj[date.vals<date+(dt/2) & date.vals>date-(dt/2)])/theta[["npop"]])*epsilon
      binom.lik[i] <- (dbinom(nLUM[i], size=nPOP[i], prob=seroP[i], log = T))
      i <- i+1
    }
  }else{
    binom.lik=0
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
  plot(date.vals, y.vals, type='l', col=4, yaxt='n', xaxt='n', xlab="", ylab="", ylim=c(0,1000))
  lines(date.vals, sapply(time.vals, function(xx){intro_f(xx, mid = theta[["zika_start_point"]], width = theta[["intro_width"]], base = theta[["intro_base"]])}),
        type='l', col=1, yaxt='n', xaxt='n', xlab="", ylab="")
  axis(side=4)
  mtext("Introductions", side=4, cex = 0.7, padj=3)
  
  # Plot seasonality
  plot(date.vals,rr_post,type='l',col=1, ylim=c(0,3), xlab="", ylab="")
  lines(date.vals,decline_ii, type='l',col=2, yaxt='n', xaxt='n')
  par(new=T)
  plot(date.vals, (R_traj/342000)+theta[["epsilon"]], type='l', col=22, yaxt='n', xaxt='n', xlab="", ylab="", ylim=c(0,1))
  points(c(as.Date("2013-11-15"), max(date.vals)), (nLUM[1:2]/nPOP[1:2]), col=4, pch=16)
  axis(side=4)
  
  return(likelihood)
}
