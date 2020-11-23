# - - - - - - - - - - - - - - - - - - - - - - - 
# Code to run one simulation of ZIKV transmission model
# Author: Alasdair Henderson
# github.com/
# - - - - - - - - - - - - - - - - - - - - - - - 

zika_single_sim <- function(transmission_rate){
  zika_data <- data.frame(date=as.Date(date.vals), Zika2016=y.vals)
  zika_data$date <- as.Date(zika_data$date)
  
  # Change parameter values -------------------------------------------------
  par(mfrow=c(3,1), mar = c(3,4,1,3))
  
  theta <- c(thetatab[1,],thetaAlltab[1,iiH,], theta_denv)
  theta_init <- theta_initAlltab[1,iiH,]
  
  #theta[["intro_mid"]] <- zika_intro
  theta[["beta_h"]] <- transmission_rate
  
  # theta[["rep"]] <- 0.00012
  # theta[["alpha"]] <- 0.78
  # theta[["epsilon"]] <- 0.062
  # theta[["intro_base"]] <- 12
  # theta[["rep"]] <- 9.8e-5
  # theta[["alpha"]] <- 0.81
  # theta[["epsilon"]] <- 0.06
  # theta[["intro_base"]] <- 4.6
  # theta[["repvol"]] <- 0.1
  
  # check that R is > 1 before introducing 
  theta_test <- as.data.frame(t(theta))
  theta_test$Inf. <- theta_test[["Inf"]]
  seasonal_min <- theta_test[["intro_mid"]]-theta_test[["intro_width"]]
  seasonal_max <- theta_test[["intro_mid"]]+theta_test[["intro_width"]]
  
  b_vary <- seasonal_f(seasonal_min:seasonal_max, 
                       0, amp=theta_test[['beta_v_amp']], mid=theta_test[['beta_v_mid']])
  R_at_intro <- calculate_r0(th_in = theta_test, sus = 1, b_vary = b_vary)$rr_out
  if(mean(R_at_intro)<1){theta[["intro_base"]] <- 0}
  
  introductions_mid <- theta[["intro_mid"]]
  introductions_width <- theta[["intro_width"]]
  introductions_base <- theta[["intro_base"]]
  ## total introductions
  total_intro <- 4 * introductions_width * introductions_base # Integral of 4*base*exp(-(time-mid)/width)/(1+exp(-(time-mid)/width))^2 over -infty/infty is 4*base*width
  
  # set initial conditions
  init1 <- c(
    s_init=theta_init[["s_init"]],e_init=theta_init[["i1_init"]],i_init=theta_init[["i1_init"]],r_init=theta_init[["r_init"]],c_init=0,
    sd_init=theta_init[["sd_init"]],ed_init=theta_init[["ed_init"]],id_init=theta_init[["id_init"]],t1d_init=theta_init[["t1d_init"]],t2d_init=theta_init[["t2d_init"]],cd_init=0)
  
  # time series, date series and case data
  # adjust when DENV and ZIKV introductions begin (depending on when model start date is): n.b. DENV3 fixed to 2013-10-27
  theta[["denv_start_point"]] <- as.Date("2013-10-29") - startdate
  
  theta[["zika_start_point"]] <- theta[["intro_mid"]]
  
  if(!is.na(theta[['epsilon']])){
    epsilon <- theta[['epsilon']]}else{
      epsilon <- 0}
  if(!is.na(theta[['alpha']])){
    alpha <- theta[['alpha']]}else{
      alpha <- 0}
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
  mu_plot <- theta[["mu"]]
  eta_plot <- theta[["eta"]]
  # Run simulation ----------------------------------------------------------
  output <- zikv_model_ode(theta, init1, time.vals.sim=time.vals)
  
  # Match compartment states at sim.vals time
  S_traj <- output[match(time.vals,output$time),"s_init"]
  I_traj <- output[match(time.vals,output$time),"i_init"]
  R_traj <- output[match(time.vals,output$time),"r_init"]
  cases1 <- output[match(time.vals,output$time),"c_init"]
  casesD <- output[match(time.vals,output$time),"cd_init"]
  casecountD <- casesD-c(0,casesD[1:(length(time.vals)-1)])
  casecount <- cases1-c(0,cases1[1:(length(time.vals)-1)])
  cases_est <- ReportC(casecount, theta[["rep"]], 1/theta[["repvol"]] )
  cases_est[cases_est==0] <- NA
  casecount[casecount<0] <- 0

  # Compute likelihood ------------------------------------------------------
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
    
  likelihood <- sum(binom.lik) + sum(log(dnbinom(y.vals,
                                                 mu=theta[["rep"]]*(casecount),
                                                 size=1/theta[["repvol"]])))
                                     
                                     # temp <- log(dnbinom(y.vals,
                                     #             mu=theta[["rep"]]*(casecount),
                                     #             size=0.001))
                                     # plot(temp, main = sum(temp))
                                     # dnbinom(1, mu = 5, size = 0.5)
                                     # dnbinom(1, mu = 0, size = 0.5)
                                     # plot(log(I_traj), ylim = c(0,10))
                                     # plot(I_traj, ylim = c(0,100), type = "l")
  # calculate R0 ------------------------------------------------------------
  tMax <- length(casecount)
  t.start = 0
  time.V = (1:tMax)*dt
  model_st <- startdate
  date0 = (model_st-date.vals[1]) %>% as.numeric() 
  
  beta_ii <- seasonal_f(time.V[1:tMax], date0, amp=theta[["beta_v_amp"]], mid=theta[["beta_v_mid"]])
    #decline_ii <- decline_f(time.V[1:tMax], mid=theta[["beta_mid"]], width=theta[["beta_grad"]], base=theta[["beta_base"]])
  decline_ii <- control_f(time.V[1:tMax], base=theta[["beta_base"]], mid=theta[["beta_mid"]], width=theta[["beta_width"]])
  b_vary = beta_ii#*decline_ii
  
  s_pick = S_traj[1:tMax]/342000
  theta[["Inf."]] <- theta[["Inf"]]
  theta[["mu"]] <- 67.25
  
  output_rr <- calculate_r0(th_in = as.data.frame(t(theta)), sus = s_pick, b_vary = b_vary)
  
  start.rr <- output_rr$rr_out[min(which(output_rr$rr_out>0))]
  output_rr$rr_out[1:(min(which(output_rr$rr_out>0))-1)] <- start.rr
  output_rr$r0_out[1:(min(which(output_rr$r0_out>0))-1)] <- start.rr
  rr_post = output_rr$rr_out
  rr_post[length(rr_post)] = rr_post[length(rr_post)-1]
  r0_post = output_rr$r0_out
  r0_post[length(r0_post)] = r0_post[length(r0_post)-1]
  decline_post= decline_ii  
  
  # Plot simulation ---------------------------------------------------------
  plot(date.vals, casecount, type='l', col=4, lwd = 2, xlab="Year", ylab="Infections (cases)",ylim = c(0, 15000),
       main=paste0(signif(sum(binom.lik),3), " Start: ", startdate+theta[["intro_mid"]],"   |   lik:", signif(likelihood,4),"  |  medR0: ", signif(median(r0_post),3), " | intro:", round(total_intro,0)))
  #points(date.vals, c(denv_data$y.vals, rep(NA, 112)), col = "gray60")
  lines(date.vals,casecountD)
  par(new=T)
  plot(date.vals, y.vals, type='h', col=4, yaxt='n', xaxt='n', xlab="", ylab="", ylim=c(0,10))
  lines(date.vals,sapply(time.vals, function(xx){intro_f(xx, mid = theta[["zika_start_point"]], width = theta[["intro_width"]], base = theta[["intro_base"]])}),
        type='l',col = 3, yaxt='n', xaxt='n', xlab="", ylab="")
  axis(side=4)
  mtext("Introductions", side=4, cex = 0.7, padj=3)
  
  # Plot seasonality
  plot(date.vals,rr_post,type='l',col=1, ylim=c(0,3), xlab="", ylab="")
  lines(date.vals,decline_ii, type='l',col=2, yaxt='n', xaxt='n')
  par(new=T)
  adjusted_pop <- sapply(time.vals, function(xx){theta[["npop"]]-(xx*mu_plot*theta[["npop"]])+(xx*eta_plot*theta[["npop"]])})
    modelled_R_traj <- R_traj/adjusted_pop
      detected_positives <- modelled_R_traj*alpha    ## identify alpha% of the actual positives
      false_positives <- (1-modelled_R_traj)*epsilon ## falsely record epsilon% of the actual negatives as positives
      
  recorded_sero <- false_positives + detected_positives
  
  plot(date.vals, recorded_sero, type='l', col=22, yaxt='n', xaxt='n', xlab="", ylab="", ylim=c(0,1))
  lines(date.vals, modelled_R_traj, type='l', col=22, lty = 2, yaxt='n', xaxt='n', xlab="", ylab="", ylim=c(0,1))
  points(seroposdates, (nLUM/nPOP), col=4, pch=1)
  axis(side=4)
  
  # plot susceptibles
  plot(date.vals, I_traj, col = 8, type = "l")
  par(new=T)
  plot(date.vals, (S_traj/342000), type='l', col=6, yaxt='n', xaxt='n', xlab="", ylab="", ylim=c(0,1.5))
  axis(side =4)
  return(likelihood)
}
