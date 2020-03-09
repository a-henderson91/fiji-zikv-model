# - - - - - - - - - - - - - - - - - - - - - - - 
# MCMC analysis of Zika transmission in Fiji  
# Author: Alasdair Henderson
# github.com/a-henderson91/fiji-zikv-model
# - - - - - - - - - - - - - - - - - - - - - - - 
library(here)
preamblecode <- "preamble_zikvfiji.R"
source(here::here("Rscripts", preamblecode))

r_functions <- list.files(here::here("Rfunctions/"), pattern="*.R")
sapply(here::here("Rfunctions", r_functions), source, .GlobalEnv)
data <- load.data.multistart(Virus = "ZIKV", startdate = start.output.date, serology.file.name = serology.excel, init.values.file.name = init.conditions.excel, add.nulls = 0) #virusTab[iiH], dataTab[iiH])

# Load and store denv2014 posteriors -----------------------------------------------
labelN=1
labelN <- 1
iiH <- 2
load_posterior_1 <- load.posteriors(load.run.name=model1_name, file.path="posterior_denv2014fit", iiH, mcmc.burn=mcmc.burn)
  list2env(load_posterior_1,globalenv())
iiH <- 1
## fixed DENV3 2014 values to simulate 2013-14 DENV-3 outbreak
load(file=here::here("data/theta_denv.RData"))

## load data
denv_data <- load.data.multistart(Virus = "DEN3", startdate = start.output.date, serology.file.name = serology.excel, init.values.file.name = init.conditions.excel, add.nulls = 0) #virusTab[iiH], dataTab[iiH])
  list2env(denv_data,globalenv())
denv.timeseries <- denv_data$y.vals
denv.dates <- denv_data$date.vals

## plot posterior of cases
tMax <- dim(c_trace_tab)[2]
btsp <- 4000
cvector <- matrix(NA,nrow=btsp,ncol=tMax)
for(ii in 1:btsp){
    pick <- sample(picks, 1)
    cvector[ii,] <- ReportC(c_trace_tab[pick,1:tMax],thetatab[pick,'rep'],thetatab[pick,'repvol'])
}
medP <- apply(cvector,2,function(x){median(x, na.rm=T)})
ciP1 <- apply(cvector,2,function(x){quantile(x,0.025, na.rm=T)})
ciP2 <- apply(cvector,2,function(x){quantile(x,0.975, na.rm=T)})

par(mfrow=c(1,1))
plot(denv.dates, denv.timeseries, ylim=c(0, max(ciP2)), bty="n", type='b', ylab="Cases", xlab="Date", pch=16, col=1)
  lines(denv.dates[1:tMax], medP, col=col1)
  polygon(c(denv.dates[1:tMax], rev(denv.dates[1:tMax])),
          c(ciP1, rev(ciP2)), lty=0, col=col1a)
mtext(LETTERS[1],side=3, adj=0, font=2)
grid(NA,NULL, lty = 1, col = colgrid) 

# calculate empirical distributions from DENV3 2014 fit -------------------
model1_pars_estimated <- apply(thetatab,2,function(x){var(x)})
priorFit <- fitdistr(thetatab$beta_base, densfun="normal")

prior_vectorcontrol <- function(x){dnorm(x, mean=priorFit$estimate[1], sd=priorFit$estimate[2])}
hist(thetatab[["beta_base"]], prob=TRUE, main="Beta_base estimate & prior", xlab="Beta_base")
curve(prior_vectorcontrol(x), col=col2, lwd=2, add=T, yaxt="n")

## save posterior estiamtes from 2014 DENV 
denv3_2014_esimates <- data.frame(NA)
denv3_2014_esimates$beta_h <- median(thetatab$beta_h, na.rm=T)
denv3_2014_esimates$beta_v <- median(thetatab$beta_v, na.rm=T)
denv3_2014_esimates$beta_v_amp <- median(thetatab$beta_v_amp, na.rm=T)
denv3_2014_esimates$beta_v_mid <- median(thetatab$beta_v_mid, na.rm=T)
denv3_2014_esimates$beta_grad <- median(thetatab$beta_grad, na.rm=T)
denv3_2014_esimates$beta_mid <- median(thetatab$beta_mid, na.rm=T)
denv3_2014_esimates$beta_base <- median(thetatab$beta_base, na.rm=T)
denv3_2014_esimates$beta_width <- median(thetatab$beta_width, na.rm=T)
denv3_2014_esimates$mosq_ratio <- median(thetatab$m, na.rm=T)

# Load ZIKV data  --------------------------------------------------------------
list2env(data,globalenv())

# Set up results objects, initial values and covariance matrices ----------
setup <- results_set_up(iiH, parameter_est_file = "parameters_est")
  list2env(setup, globalenv())

# Set up Priors  ----------------------------------------------------------
## Estimate prior distribution for ZIKV intro time from BEAST output
source(here::here("Rscripts/load_tmrca_calc_prior.R"))
intro_prior_mu;intro_prior_sigma
source(here::here("Rscripts/prior_distributions.R"))

# Set initial values to posterior means form denvlike fit -----------------
thetaAlltab[1,iiH,'beta_v_amp'] <- denv3_2014_esimates$beta_v_amp ## seasonality parameters estimated during DENV3-2014 fitting
thetaAlltab[1,iiH,'beta_v_mid'] <- denv3_2014_esimates$beta_v_mid ## seasonality parameters estimated during DENV3-2014 fitting
thetaAlltab[1,iiH,'beta_grad'] <- denv3_2014_esimates$beta_grad
thetaAlltab[1,iiH,'beta_mid'] <- denv3_2014_esimates$beta_mid
thetaAlltab[1,iiH,'beta_width'] <- denv3_2014_esimates$beta_width
thetaAlltab[1,iiH,'m'] <- denv3_2014_esimates$mosq_ratio

if(include.2014.control == T){
  thetaAlltab[1,iiH,'beta_base'] <- denv3_2014_esimates$beta_base
}else{
  thetaAlltab[1,iiH,'beta_base'] <- 0
}

## Plot "control function" with initial starting values
control_plot_vals <- as.data.frame(t(thetaAlltab[1,iiH,]))
control_plot <- sapply(time.vals, function(xx){control_f(xx, 
                                         base=control_plot_vals$beta_base,
                                         grad=control_plot_vals$beta_grad,
                                         mid=control_plot_vals$beta_mid,
                                         mid2=control_plot_vals$beta_mid+control_plot_vals$beta_width,
                                         width=control_plot_vals$beta_width
                                         )})
plot(date.vals[30:100], control_plot[30:100], type="l", ylim=c(0,1), bty="n", ylab="Relative R0", xlab="Date", xaxt="n")
axis.Date(side=1, at=seq.Date(date.vals[30], date.vals[100], by = "1 months"), "months", format = "%b%y")

dev.copy(pdf, here::here("output/fig_supp_controlFn.pdf"), 6,4)
  dev.off()
  
# run a single simulation with initial values -----------------------------
source(here::here("Rscripts/zika-single-simulation.R"))
initial_beta_h <- thetaAlltab[1,iiH,][["beta_h"]]
zika_single_sim(initial_beta_h)
if(zika_single_sim(initial_beta_h)==-Inf){
  t_rate <- seq(0.3,0.4,0.01)
  lik_search <- sapply(seq(0.3,0.4,0.01), function(xx){zika_single_sim(xx)})
  max_beta_h <- t_rate[which(lik_search==max(lik_search))]
  thetaAlltab[1,iiH,"beta_h"] <- max_beta_h
}

## Plot "introduction function" with initial starting values
intro_plot_vals <- as.data.frame(t(thetaAlltab[1,iiH,]))
intro_plot <- sapply(time.vals, function(xx){1-intro_f(xx, 
                                                         mid=control_plot_vals$beta_mid,
                                                         width=control_plot_vals$beta_width,
                                                         base=control_plot_vals$intro_base
)})
plot(date.vals[20:90], intro_plot[20:90], type="l", ylim=c(0,1), bty="n", ylab="ZIKV introductions", xlab="Date", xaxt="n")
axis.Date(side=1, at=seq.Date(date.vals[20], date.vals[90], by = "1 months"), "months", format = "%b%y")

dev.copy(pdf, here::here("output/fig_supp_introFn.pdf"), 6,4)
  dev.off()
  
# Print starting conditions for model run --------------------------
print("Theta initial starting conditions")
thetaAlltab[1,iiH,]
diag(cov_matrix_thetaAll)

# Run MCMC loop and save results ------------------------------------------
(st.time <- Sys.time())

m = 1; iiM=1
foreach(iiM=multichain) %dopar% {  # Loop over regions with parallel MCMC chains

adapt_size_start <- round(0.1*MCMC.runs)

for (m in 1:MCMC.runs){
  # Scale COV matrices for resampling using error term epsilon
  if(m==1){
    epsilon0 = 0.001
    accept_rate = NA
    cov_matrix_theta=epsilon0*cov_matrix_theta0
    cov_matrix_thetaA=epsilon0*cov_matrix_thetaAll
    cov_matrix_theta_init=epsilon0*cov_matrix_theta_initAll
    # set current values for vectors to be updated in MH algroithm 
    thetatab_current = thetatab[m,]
    thetaAlltab_current = thetaAlltab[m,,]
    theta_initAlltab_current = theta_initAlltab[m,,]
    c_trace_tab_current = c_trace_tab[m,,]
    cd_trace_tab_current = cd_trace_tab[m,,]
    s_trace_tab_current = s_trace_tab[m,,]
    r_trace_tab_current = r_trace_tab[m,,]
    x_trace_tab_current =  x_trace_tab[m,,]
    sim_liktab_current = sim_liktab[m]
    prior_current = prior[m]
    covmat.empirical <- cov_matrix_thetaA
    theta.mean <- thetaAlltab_current[1,]
    adapting.shape <- 0
    # initialise counter for storing results (m/thining parameter)
    j=1
  }else{
    scaling.multiplier <- exp((1-1e-7)^(m-adapt_size_start)*(accept_rate-0.234))
    epsilon0 <- epsilon0 * scaling.multiplier
    epsilon0 <- min(epsilon0, 0.5)
    epsilon0 <- max(epsilon0, 1e-10)
    
    cov_matrix_theta=epsilon0*cov_matrix_theta0
    cov_matrix_thetaA=epsilon0*cov_matrix_thetaAll
    cov_matrix_theta_init=epsilon0*cov_matrix_theta_initAll
  }
  
  ## Resample global theta every 2nd step
  if(m %% 2==0){
    output_theta <- SampleTheta(global=1,thetatab_current,theta_initAlltab_current[iiH,],cov_matrix_theta,cov_matrix_theta_init)
    theta_star <- output_theta$thetaS
  }else{
    theta_star <- thetatab_current
  }
  
  ## Resample local parameters every step
  prior.star <- 1
  prior.current <- 1
  sim_marg_lik_star <- 0
  thetaAllstar <- 0*thetaAlltab_current
  theta_initAllstar <- 0*theta_initAlltab_current
  cTraceStar <- 0*c_trace_tab_current
  cdTraceStar <- 0*cd_trace_tab_current
  sTraceStar <- 0*s_trace_tab_current
  rTraceStar <- 0*r_trace_tab_current
  xTraceStar <- 0*x_trace_tab_current
  for(kk in 1){ #itertab){ 
    iiH <- kk
    data <- load.data.multistart(Virus = virusTab[iiH], startdate = start.output.date, serology.file.name = serology.excel, init.values.file.name = init.conditions.excel, add.nulls = 0) #virusTab[iiH], dataTab[iiH])
      list2env(data,globalenv())
      
    if(m==1){ # Don't resample on 1st step - check the zeroes!
      output_H <- SampleTheta(global=0,thetaAlltab_current[iiH,],theta_initAlltab_current[iiH,],0*cov_matrix_thetaA,0*cov_matrix_theta_init)
    }else{
      output_H <- SampleTheta(global=0,thetaAlltab_current[iiH,],theta_initAlltab_current[iiH,],cov_matrix_thetaA,cov_matrix_theta_init)
    } 
    thetaA_star <- output_H$thetaS
    theta_init_star <- output_H$theta_initS
    
    # Run model simulation
    output1 <- Deterministic_modelR_final_DENVimmmunity(theta=c(theta_star,thetaA_star,theta_denv), theta_init_star, locationI=locationtab[iiH], seroposdates=seroposdates, include.count=include.count)
      (sim_marg_lik_star=sim_marg_lik_star + output1$lik)
    
    #Store vales
    thetaAllstar[iiH,] <- thetaA_star
    theta_initAllstar[iiH,] <- theta_init_star
    
    # store results between prespecified dates
    start.of.output1 <- 1
    length.of.output1 <- length(output1$C_trace) - max(0,length(output1$C_trace) - length(cTraceStar[iiH,]))
    extra.vals <- max(0,length(cTraceStar[iiH,]) - length(output1$C_trace))
    extra.zero <- rep(0,extra.vals)
    cTraceStar[iiH,]=c(extra.zero,output1$C_trace[(start.of.output1:length.of.output1)])
    cdTraceStar[iiH,]=c(extra.zero,output1$CD_trace[(start.of.output1:length.of.output1)])
    sTraceStar[iiH,]=c(extra.zero,output1$S_trace[(start.of.output1:length.of.output1)])
    rTraceStar[iiH,]=c(extra.zero,output1$R_trace[(start.of.output1:length.of.output1)])
    xTraceStar[iiH,]=c(extra.zero,output1$X_trace[(start.of.output1:length.of.output1)])
    
    # Calculate prior density for current and proposed theta set
    prior.theta <- ComputePrior(iiH, c(thetatab_current,thetaAlltab_current[iiH,]), c(theta_star,thetaA_star), covartheta = cov_matrix_thetaA)
    prior.star <- (prior.theta$prior.star*prior.star)
    prior.current <- (prior.theta$prior*prior.current)
  } # end loop over regions
  
  # Calculate probability function - MH algorithm
    estimated_params <- colnames(cov_matrix_thetaA)[diag(cov_matrix_thetaA)!=0]
    q_theta_given_theta_star = sum(log(thetaAllstar[iiH, estimated_params]))
    q_theta_star_given_theta = sum(log(thetaAlltab_current[1, estimated_params]))
  
  val <- (prior.star/prior.current)*exp((sim_marg_lik_star-sim_liktab_current) +  (q_theta_given_theta_star - q_theta_star_given_theta)) 
  
  if(is.na(val)){
    output_prob=0}else if(is.nan(val)){
      output_prob=0}else if(is.null(val)){
        output_prob=0}else if(length(val)==0){
          output_prob=0}else{
            output_prob = min(val, 1)}
  
  # Update parameter values every k step
  MH_random_unif <- runif(1,0,1)
  
  if(m %% thinning.parameter == 0){
    if(MH_random_unif < output_prob){
      thetatab[j+1,] = theta_star
      thetaAlltab[j+1,,] = thetaAllstar
      theta_initAlltab[j+1,,] = theta_initAllstar
      c_trace_tab[j+1,,]=cTraceStar
      cd_trace_tab[j+1,,]=cdTraceStar
      s_trace_tab[j+1,,]=sTraceStar
      r_trace_tab[j+1,,]=rTraceStar
      x_trace_tab[j+1,,]=xTraceStar
      sim_liktab[j+1] = sim_marg_lik_star
      accepttab[j]=1
      prior[j+1] = prior.star
    }else{
      thetatab[j+1,] = thetatab[j,]
      thetaAlltab[j+1,,] = thetaAlltab[j,,]
      theta_initAlltab[j+1,,] = theta_initAlltab[j,,]
      c_trace_tab[j+1,,]=c_trace_tab[j,,]
      cd_trace_tab[j+1,,]=cd_trace_tab[j,,]
      s_trace_tab[j+1,,]=s_trace_tab[j,,]
      r_trace_tab[j+1,,]=r_trace_tab[j,,]
      x_trace_tab[j+1,,]=x_trace_tab[j,,]
      sim_liktab[j+1] = sim_liktab[j]
      accepttab[j]=0
      prior[j+1] = prior[j] 
    }
    accept_rate=sum(accepttab[1:j])/j
    j <- j+1
  }
  
  # Update current values of parameter if{MH_algorithm_val > runif(1,0,1)}
  if(MH_random_unif < output_prob){
    thetatab_current = theta_star
    thetaAlltab_current = thetaAllstar
    theta_initAlltab_current = theta_initAllstar
    c_trace_tab_current = cTraceStar
    cd_trace_tab_current = cdTraceStar
    s_trace_tab_current = sTraceStar
    r_trace_tab_current = rTraceStar
    x_trace_tab_current = xTraceStar
    sim_liktab_current = sim_marg_lik_star
    prior_current = prior.star
  }
  
  if(m<adapt_size_start){
    accept_rate <- 0.234
  }
  
  if(m %% min(MCMC.runs,1000)==0){
    print(c(iiM,"m"=m,  
            "acc"=signif(accept_rate,3), 
            "lik"=signif(sim_liktab_current,3),
            signif(thetaAlltab_current[1,'beta_h'],3),
            thetaAlltab_current[1,'intro_mid'],
            thetaAlltab_current[1,'rep']),
            digits = 2)
  save(sim_liktab,
            prior,
            accepttab,
            c_trace_tab,
            cd_trace_tab,
            s_trace_tab,
            r_trace_tab,
            x_trace_tab,
            thetatab,
            thetaAlltab,
            theta_initAlltab,
            file=paste("posterior_output/outputR_",iiM,"_",run.name,".RData",sep=""))
    }  
  } # End MCMC loop
}
endtime <- Sys.time()
(time.take = endtime-st.time)


## plot posterior of cases
source(here::here("Rscripts", preamblecode))
data <- load.data.multistart(Virus = "ZIKV", startdate = start.output.date, serology.file.name = serology.excel, init.values.file.name = init.conditions.excel, add.nulls = 0) #virusTab[iiH], dataTab[iiH])
  list2env(data,globalenv())

load_posterior_1 <- load.posteriors(load.run.name=run.name, file.path="posterior_output", iiH, mcmc.burn=mcmc.burn)
  list2env(load_posterior_1,globalenv())

tMax <- dim(c_trace_tab)[2]
btsp <- 4000
cvector <- matrix(NA,nrow=btsp,ncol=tMax)
for(ii in 1:btsp){
  pick <- sample(1:MCMC.runs, 1)
  cvector[ii,] <- ReportC(c_trace_tab[pick, 1:tMax],thetatab[pick,'rep'],thetatab[pick,'repvol'])
}
medP <- apply(cvector,2,function(x){median(x, na.rm=T)})
ciP1 <- apply(cvector,2,function(x){quantile(x,0.025, na.rm=T)})
ciP2 <- apply(cvector,2,function(x){quantile(x,0.975, na.rm=T)})

par(mfrow=c(1,1))
plot(data$date.vals, data$y.vals, ylim=c(0, max(ciP2)), bty="n", type='h', ylab="Cases", xlab="Date", pch=16, col=1)
lines(data$date.vals[1:tMax], medP, col=col1)
polygon(c(data$date.vals[1:tMax], rev(data$date.vals[1:tMax])),
        c(ciP1, rev(ciP2)), lty=0, col=col1a)
mtext(LETTERS[1],side=3, adj=0, font=2)
grid(NA,NULL, lty = 1, col = colgrid) 
