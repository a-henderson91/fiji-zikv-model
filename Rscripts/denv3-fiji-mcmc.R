# - - - - - - - - - - - - - - - - - - - - - - - 
# MCMC analysis of 2013-14 DENV-3 transmission in Fiji  
# Author: Alasdair Henderson
# github.com/a-henderson91/fiji-zikv-model
# - - - - - - - - - - - - - - - - - - - - - - - 
library(here)
preamblecode <- "preamble_denv3.R"
source(here::here("Rscripts", preamblecode))

r_functions <- list.files(here::here("Rfunctions/"), pattern="*.R")
sapply(here::here("Rfunctions", r_functions), source, .GlobalEnv)

# load data and setup with fixed DENV intro date Oct 27th 2013
prior_vectorcontrol <- function(x){
  dnorm(x, mean = 0.57, sd = 0.15) ## based on Kucharski et al. 2018 Elife
}
denv_data <- load.data.multistart(Virus = "DEN3", startdate = start.output.date, serology.file.name = serology.excel, init.values.file.name = init.conditions.excel, add.nulls = 0)
  list2env(denv_data,globalenv())
denv.timeseries <- denv_data$y.vals
denv.dates <- denv_data$date.vals

load(file=here::here("data/theta_denv.RData"))
theta_denv

# Set up results objects, initial values and covariance matrices ----------
setup <- results_set_up(iiH, "parameters_est_denv3")
  list2env(setup, globalenv())
theta_initAlltab[1,,]
thetaAlltab[1,,]


# Set up Priors  ----------------------------------------------------------
source(here::here("Rscripts/prior_distributions.R"))

priorIntro <- function(x){dunif(x, min=0, max=307)} ## DENV import mid point is 2013-11-07 (first case reported) at latest
priorInitInf <- function(x){dunif(x, min=0, max=1000)} ## DENV import max is 1000 in one day
priorEpsilon <- function(x){dtruncnorm(x, a=0, b=1, mean=0, sd=0.15)} ## 100% sensitivity and specificity in control panel 
priorAlpha <- function(x){dtruncnorm(x, a=0, b=1, mean=1, sd=0.15)}

# load climate data -------------------------------------------------------
weather.data <- read.csv(here::here("data/suva-temp.csv"), stringsAsFactors = T) # Load DLI dengue data
weather.data$lsd <- as.Date(weather.data$lsd)

first.temp.date <- min(as.Date(date.vals))
data.series.max <- weather.data$max_air_temp[as.Date(weather.data$lsd)>=first.temp.date & as.Date(weather.data$lsd)<=max(as.Date(date.vals))]
data.series.min <- weather.data$min_air_temp[as.Date(weather.data$lsd)>=first.temp.date & as.Date(weather.data$lsd)<=max(as.Date(date.vals))]
data.series <- ((data.series.max-data.series.min)/2)+data.series.min
date.series <- weather.data$lsd[as.Date(weather.data$lsd)>=first.temp.date & as.Date(weather.data$lsd)<=max(as.Date(date.vals))]

#betas <- weather.fit(data=data.series, iter=min(MCMC.runs, 5e4), tuning=10)
#  save(betas,file=paste("posterior_output_betas/outputBetas.RData",sep=""))
  load(file=paste("posterior_output_betas/outputBetas.RData",sep=""))

  plot(betas$lik[2,],betas$mid,type='p')
  max(betas$lik[2,],na.rm=T)
  
# Find seasonality parameters ---------------------------------------------
mordecai_tempR0 <- read.csv(here::here("data/outputTemp_R0.csv"))
btstrp <- 4000
sin_sum <- NULL
picks <- 1:length(betas$amp)
b <- mean(data.series, na.rm = T)
for(ii in 1:btstrp){
  b_ii <- sample(picks,1)
  sincurve <- b*seasonal_f(time = time.vals, date0 = 0, amp = betas$amp[b_ii], mid = betas$mid[b_ii])
  sin_sum <- rbind(sin_sum, sincurve)
}
c.nume<-function(x){
  bp1=c(median(x),quantile(x,0.25),quantile(x,0.75))
  as.numeric(bp1)}
sinwave_summ <- apply(sin_sum,2,c.nume) 
plot(date.series, data.series, col="gray75", type='l')
lines(date.vals, sinwave_summ[1,], type='l')
peakBeta <- max(sinwave_summ[1,])
lowBeta <- min(sinwave_summ[1,])
peakBetaDate <- date.vals[sinwave_summ[1,]==peakBeta]
lowBetaDate <- date.vals[sinwave_summ[1,]==lowBeta]
peakBetaTemp <- round(peakBeta,1) #round(sum(weather.data[weather.data$lsd==peakBetaDate,3:4])/2,1)
lowBetaTemp <- round(lowBeta,1) #round(sum(weather.data[weather.data$lsd==lowBetaDate,3:4])/2,1)
## find in Mordecai data
(peakR0 <- mordecai_tempR0[round(mordecai_tempR0$aegy.temps.DTR8,1)==peakBetaTemp,2])
(lowR0 <- mordecai_tempR0[round(mordecai_tempR0$aegy.temps.DTR8,1)==lowBetaTemp,2])
rel_lowR0 <- lowR0/peakR0 
#print(rel_lowR0)
mordecai_amp <- (1-rel_lowR0)/(1+rel_lowR0)
print(c("mordecai_amp",mordecai_amp))

# Show Fiji temperature on Mordecai data ----------------------------------
med_minT <- median(weather.data$min_air_temp, na.rm=T)
lq_minT <- quantile(weather.data$min_air_temp, 0.25, na.rm=T)
uq_minT <- quantile(weather.data$min_air_temp, 0.75, na.rm=T)
med_maxT <- median(weather.data$max_air_temp, na.rm=T)
lq_maxT <- quantile(weather.data$max_air_temp, 0.25, na.rm=T)
uq_maxT <- quantile(weather.data$max_air_temp, 0.75, na.rm=T)
med_avgT <- median(data.series, na.rm=T)
lq_avgT <- quantile(data.series, 0.25, na.rm=T)
uq_avgT <- quantile(data.series, 0.75, na.rm=T)

## plot
par(mfrow = c(1,2))
plot(date.series, data.series, type='l', col = col1, ylab = "Temperature (c)", xlab = "Year", bty = "n")
grid(NA,NULL, lty = 1, col = colgrid) 
mtext("A",side=3, adj=0, font=2)
plot(mordecai_tempR0$aegy.temps.DTR8, mordecai_tempR0$R0.rel, type="l", bty="n", 
     xlab="Temperature (C)", ylab="Relative R0")
lines(c(med_avgT,med_avgT), c(0,1), col=col1)
  polygon(c(lq_avgT, uq_avgT, uq_avgT, lq_avgT), c(0,0,1,1), lty=0, col=col1a)
mtext("B",side=3, adj=0, font=2)
dev.copy(pdf, here::here("output/EDITED_fig_supp_mordecaiFijitemp.pdf"), 6, 4)
  dev.off()

# Check seasonality fit worked --------------------------------------------
beta=mean(data.series, na.rm = T)
(amp=median(betas$amp))
(amp=mordecai_amp)
(mid=median(betas$mid)) #(mid=0.00027)
j=1; sincurve=NULL; date0=0
for(i in 1:length(date.series)){
  sincurve[j] <- beta*(1 + amp * sin(((i/365.25) - mid)*2*pi))
  j=j+1
}
plot(date.series, data.series, type='l', ylim=c(min(min(data.series, na.rm=T), min(sincurve)), max(max(data.series, na.rm=T), max(sincurve))))
lines(date.series, sincurve, col=2, type='l',yaxt="n")
max_curve <- date.series[sincurve==max(sincurve)]
print(paste0("max curve = ", max_curve))
if(months(max_curve)!="February"){
  error("Peak month of seeasonal wave is not in Feb. Check why")
}

burnin=0.2
if(seasonal.transmission == T){
  thetaAlltab[1,,'beta_v_amp'] <- mordecai_amp
  thetaAlltab[1,,'beta_v_mid'] <- median(betas$mid[round(burnin*length(betas$mid)):length(betas$mid)])
}else{
  thetaAlltab[1,,'beta_v_amp'] <- 0
}  

# set weakly informative priors for main model fit 
priorBeta_amp <- function(x){dunif(x,min=0, max=1)} 
priorBeta_mid <- function(x){dunif(x,min=0, max=1)} 

# Print starting conditions for model run --------------------------
print("Theta initial starting conditions")
thetaAlltab[1,,]
diag(cov_matrix_thetaAll)

iiH <- 2
source(here::here("Rscripts/den3-single-simulation.R"))
initial_beta_h <- thetaAlltab[1,iiH,][["beta_h"]]
denv3_single_sim(initial_beta_h)
#if(denv3_single_sim(initial_beta_h) == -Inf){
  t_rate <- seq(0.3, 0.6, 0.01)
  lik_search <- sapply(t_rate, function(xx){denv3_single_sim(xx)})
  max_beta_h <- t_rate[which(lik_search==max(lik_search))]
  thetaAlltab[1,iiH,"beta_h"] <- max_beta_h[1]
denv3_single_sim(thetaAlltab[1,iiH,"beta_h"])
#}
thetaAlltab[1,,]
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
      cov_matrix_theta <- epsilon0*cov_matrix_theta0
      cov_matrix_thetaA <- epsilon0*cov_matrix_thetaAll
      cov_matrix_theta_init <- epsilon0*cov_matrix_theta_initAll
      # set current values for vectors to be updated in MH algroithm 
      thetatab_current = thetatab[m,]
      thetaAlltab_current = thetaAlltab[m,,]
      theta_initAlltab_current = theta_initAlltab[m,,]
      c_trace_tab_current = c_trace_tab[m,,]
      cd_trace_tab_current = cd_trace_tab[m,,]
      s_trace_tab_current = s_trace_tab[m,,]
      r_trace_tab_current = r_trace_tab[m,,]
      sim_liktab_current = sim_liktab[m]
      prior_current = prior[m]
      covmat.empirical <- cov_matrix_thetaA
      theta.mean <- thetaAlltab_current[1,]
      # initialise counter for storing results (m/thining parameter)
      j=1
    }else{
      scaling.multiplier <- exp((0.9)^(m-adapt_size_start)*(accept_rate-0.234))
      epsilon0 <- epsilon0 * scaling.multiplier
      #epsilon0 <- min(epsilon0, 0.5)
      #epsilon0 <- max(epsilon0, 1e-15)
      
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
    prior.star=1
    prior.current=1
    sim_marg_lik_star=0
    thetaAllstar=0*thetaAlltab_current
    theta_initAllstar=0*theta_initAlltab_current
    cTraceStar=0*c_trace_tab_current
    cdTraceStar=0*cd_trace_tab_current
    sTraceStar=0*s_trace_tab_current
    rTraceStar=0*r_trace_tab_current
    for(kk in 2){ ## settng kk==2 means running the model for DENV3 not ZIKV
      iiH <- kk
      data <- load.data.multistart(Virus = "DEN3", startdate = start.output.date, serology.file.name = serology.excel, init.values.file.name = init.conditions.excel, add.nulls = 0) #virusTab[iiH], dataTab[iiH])
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
          (sim_marg_lik_star <- sim_marg_lik_star + output1$lik)
        
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
        
        # Calculate prior density for current and proposed theta set
        prior.theta <- ComputePrior(iiH, c(thetatab_current,thetaAlltab_current[iiH,]), c(theta_star,thetaA_star), covartheta = cov_matrix_thetaA)
        prior.star <- (prior.theta$prior.star*prior.star)
        prior.current <- (prior.theta$prior*prior.current)
    } # end loop over regions
    
    # Calculate probability function - MH algorithm
    estimated_params <- colnames(cov_matrix_thetaA)[diag(cov_matrix_thetaA)!=0]
    q_theta_given_theta_star <- sum(log(thetaAllstar[2, estimated_params]))
    q_theta_star_given_theta <- sum(log(thetaAlltab_current[2, estimated_params]))
    
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
        thetatab[j+1,] <- theta_star
        thetaAlltab[j+1,,] <- thetaAllstar
        theta_initAlltab[j+1,,] <- theta_initAllstar
        c_trace_tab[j+1,,] <- cTraceStar
        cd_trace_tab[j+1,,] <- cdTraceStar
        s_trace_tab[j+1,,] <- sTraceStar
        r_trace_tab[j+1,,] <- rTraceStar
        sim_liktab[j+1] <- sim_marg_lik_star
        accepttab[j] <- 1
        prior[j+1] <- prior.star
      }else{
        thetatab[j+1,] = thetatab[j,]
        thetaAlltab[j+1,,] = thetaAlltab[j,,]
        theta_initAlltab[j+1,,] = theta_initAlltab[j,,]
        c_trace_tab[j+1,,]=c_trace_tab[j,,]
        cd_trace_tab[j+1,,]=cd_trace_tab[j,,]
        s_trace_tab[j+1,,]=s_trace_tab[j,,]
        r_trace_tab[j+1,,]=r_trace_tab[j,,]
        sim_liktab[j+1] = sim_liktab[j]
        accepttab[j]=0
        prior[j+1] = prior[j] 
      }
      accept_rate <- sum(accepttab[1:j])/j
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
      sim_liktab_current = sim_marg_lik_star
      prior_current = prior.star
    }
    
    if(m<adapt_size_start){
      accept_rate <- 0.234
    }
    
    if(m %% min(MCMC.runs,1000)==0){
      print(c(iiM,"m" = m,  
              "acc" = signif(accept_rate,3), 
              "lik" = signif(sim_liktab_current,3),
              signif(thetaAlltab_current[2,'beta_h'],3),
              thetaAlltab_current[2,'beta_base'],
              thetaAlltab_current[2,'intro_base']),
              digits = 2)
      save(sim_liktab,
           prior,
           accepttab,
           c_trace_tab,
           cd_trace_tab,
           s_trace_tab,
           r_trace_tab,
           thetatab,
           thetaAlltab,
           theta_initAlltab,
           file=paste("posterior_denv2014fit/outputR_",iiM,"_",run.name,".RData",sep=""))
    }
  } # End MCMC loop
}

endtime <- Sys.time()
(time.take = endtime-st.time)