#' Deterministic SEIR-SEI function with no age structure
#'
#' This function solves an SEIR-SEI vector borne disease model.
#' @param theta Vector of parameters for model
#' @param theta_init Vector of initial conditions for model
#' @param locationI Name of location of data
#' @param seroposdates Vector with dates of seroprevalence surveys
#' @param episeason Vector with start and end point of epidemic case data
#' @param include.count True or False - whether to include count data in likelihood. Defaults to True
#' @keywords deterministic
#' @export
#theta=c(theta_star,thetaA_star,theta_denv); theta_init=theta_init_star; locationI=locationtab[iiH]; seroposdates=seroposdates; include.count=include.count
Deterministic_modelR_final_DENVimmmunity <- function(theta, theta_init, locationI, seroposdates, episeason, include.count=T){
    theta[["denv_start_point"]] <- as.Date("2013-10-27")-startdate ## 
    theta[["zika_start_point"]] <- theta[["intro_mid"]] ## these are poorly named. "zika_start_point" just refers to the main disease of interest. And 'denv_start_point' is the fixed background DENV3

    # These values tell how to match states of compartment with data points
    sim.vals <- seq(0,max(time.vals)-min(time.vals),7) + 7
    time.vals.sim <- seq(0,max(sim.vals),dt)

    # set initial conditions
    init1=c(
      s_init=theta_init[["s_init"]],e_init=theta_init[["i1_init"]],i_init=theta_init[["i1_init"]],r_init=theta_init[["r_init"]],c_init=0,
      sd_init=theta_init[["sd_init"]],ed_init=theta_init[["ed_init"]],id_init=theta_init[["id_init"]],t1d_init=theta_init[["t1d_init"]],t2d_init=theta_init[["t2d_init"]],cd_init=0,
      sm_init=theta_init[["sm_init"]],em_init=theta_init[["em_init"]],im_init=theta_init[["im_init"]])

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

    # Output simulation data
    output <- simulate_deterministic_noage_DENVimm(theta, init1, time.vals.sim)
    
    # Match compartment states at sim.vals time
    S_traj <- output[match(time.vals.sim,output$time),"s_init"]
    X_traj <- output[match(time.vals.sim,output$time),"sm_init"]
    R_traj <- output[match(time.vals.sim,output$time),"r_init"]
    I_traj <- output[match(time.vals.sim,output$time),"i_init"]
    cases1 <- output[match(time.vals.sim,output$time),"c_init"]
    casesD <- output[match(time.vals.sim,output$time),"cd_init"]
    SD <- output[match(time.vals.sim,output$time),"sd_init"]
    ED <- output[match(time.vals.sim,output$time),"ed_init"]
    ID <- output[match(time.vals.sim,output$time),"id_init"]
    RD <- output[match(time.vals.sim,output$time),"rd_init"]
    casecountD <- casesD-c(0,casesD[1:(length(time.vals.sim)-1)])
    casecount <- cases1-c(0,cases1[1:(length(time.vals.sim)-1)])
    casecount[casecount<0] <- 0
      #plot(date.vals, casecount[1:length(date.vals)])
      #plot(date.vals[1:length(casecount)],ReportC(cases = casecount,rep = theta['rep'], repvol = theta['repvol']),type='l', col=4)
      #points(date.vals,y.vals,type='l',col=2)
      #######
      #plot(date.vals[1 :length(casecount)],casecount,type='l',col=4)
      #abline(v=startdate+theta[["zika_start_point"]])
      #par(new=T)
      #plot(date.vals, sapply(time.vals, function(xx){control_f(xx, base=theta[["beta_base"]], grad=0.25, mid=theta[["beta_mid"]], mid2=theta[["beta_mid"]]+30, width=theta[["beta_grad"]])}), type='l', ylab="", yaxt="n", ylim=c(0,1))
      #par(new=T)
      #plot(date.vals[1:length(casecount)],R_traj/theta[["npop"]],type='l',col=4,yaxt='n',xaxt='n',ylim=c(0,1))
      #axis(side=4)

    # Calculate seropositivity at pre-specified dates and corresponding likelihood
    i=1; seroP=NULL; binom.lik=NULL
    sero.years <- format(as.Date(seroposdates, format="%d/%m/%Y"),"%Y")
    sero.y <- substr(sero.years,3,4)
    lum.y <- c("13","15","17")
    if(include.sero.likelihood==T){
    for(date in seroposdates){
      if(date < min(date.vals)){ # if seroprevalence date is before Zika timeline { seroprevalence = Luminex data + epislon}
        seroP[i] <- epsilon
        binom.lik[i] <- (dbinom(nLUM[lum.y==sero.y[i]], size=nPOP[lum.y==sero.y[i]], prob=seroP[i], log = T))
      }else{ # else { seroprevalence = model predicted recovered as a proportion of pop + epsilon }
          seroP[i] <-  (min(R_traj[date.vals<date+3.5 & date.vals>date-3.5])/theta[["npop"]]) +
                        (1 - min(R_traj[date.vals<date+3.5 & date.vals>date-3.5])/theta[["npop"]])*epsilon
          binom.lik[i] <- (dbinom(nLUM[lum.y==sero.y[i]], size=nPOP[lum.y==sero.y[i]], prob=seroP[i], log = T))
          }
        i <- i+1
      }
    }else{
      binom.lik=0
    }
    date <- seroposdates[1]
    ln.denv <- length(denv.timeseries)
    ln.full <- length(y.vals)
    first.zikv <- min(which(y.vals>0))
    likelihood <- sum(binom.lik) + sum(log(dnbinom(y.vals,
                                                    mu=theta[["rep"]]*(casecount),
                                                   size=1/theta[["repvol"]])))
    likelihood=max(-1e10, likelihood)
      if(is.null(likelihood)){likelihood=-1e10}
      if(is.nan(likelihood)){likelihood=-1e10}
      if(is.na(likelihood)){likelihood=-1e10}
      if(length(likelihood)==0){likelihood=-1e10}
      if(likelihood == -Inf){likelihood=-1e10}

    # Return results
    output1=list(C_trace=casecount,CD_trace=casecountD,I_trace=I_traj,
                 S_trace=S_traj,R_trace=R_traj,X_trace=X_traj,
                 lik=likelihood, newDates=date.vals)
}
