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
data <- load.data.multistart.month(Virus = "ZIKV", startdate = start.output.date, serology.file.name = serology.excel, init.values.file.name = init.conditions.excel, add.nulls = 0) 

# Load and store denv2014 posteriors -----------------------------------------------
labelN <- 1
iiH <- 2
load_posterior_1 <- load.posteriors(load.run.name = model1_name, file.path="posterior_denv2014fit", iiH, mcmc.burn=mcmc.burn)
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
tMax <- length(load_posterior_1$c_trace_tab[1,])
btsp <- 4000
cvector <- matrix(NA,nrow=btsp,ncol=tMax)
for(ii in 1:btsp){
    pick <- sample(picks, 1)
    cvector[ii,] <- ReportC(load_posterior_1$c_trace_tab[pick,1:tMax],thetatab[pick,'rep'],thetatab[pick,'repvol'])
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
priorFit <- fitdistr(thetatab$beta_base, densfun="normal")

prior_vectorcontrol <- function(x){dnorm(x, mean=priorFit$estimate[1], sd=priorFit$estimate[2])}
  hist(thetatab[["beta_base"]], prob=TRUE, main="Beta_base estimate & prior", xlab="Beta_base")
  curve(prior_vectorcontrol(x), col=col2, lwd=2, add=T, yaxt="n")

## save posterior estiamtes from 2014 DENV 
denv3_2014_esimates <- data.frame(NA)
denv3_2014_esimates$beta_h <- median(thetatab$beta_h, na.rm=T)
denv3_2014_esimates$beta_v_amp <- median(thetatab$beta_v_amp, na.rm=T)
denv3_2014_esimates$beta_v_mid <- median(thetatab$beta_v_mid, na.rm=T)
denv3_2014_esimates$beta_grad <- median(thetatab$beta_grad, na.rm=T)
denv3_2014_esimates$beta_mid <- median(thetatab$beta_mid, na.rm=T)
denv3_2014_esimates$beta_base <- median(thetatab$beta_base, na.rm=T)
denv3_2014_esimates$beta_width <- median(thetatab$beta_width, na.rm=T)

# Load ZIKV data  --------------------------------------------------------------
list2env(data,globalenv())

# Set up results objects, initial values and covariance matrices ----------
setup <- results_set_up_monthly(iiH, parameter_est_file = "parameters_est")
  list2env(setup, globalenv())
  setup$thetaAlltab[1,,]
  setup$thetatab[1,]
  setup$theta_initAlltab[1,1,]
  setup$cov_matrix_theta0
  
# Set up Priors  ----------------------------------------------------------
## Estimate prior distribution for ZIKV intro time from BEAST output
source(here::here("Rscripts/load_tmrca_calc_prior.R"))
  intro_prior_mu; intro_prior_sigma
source(here::here("Rscripts/prior_distributions.R"))

if(limit.to.2013 == T){
  priorIntro <- function(x){dunif(x, min = 1, max = 365)}
  thetaAlltab[1,iiH,'intro_mid']  <- 120
  thetaAlltab[1,iiH,'intro_base']  <- 10
}

# Set initial values to posterior means form denvlike fit -----------------
thetaAlltab[1,iiH,'beta_v_amp'] <- denv3_2014_esimates$beta_v_amp ## seasonality parameters estimated during DENV3-2014 fitting
thetaAlltab[1,iiH,'beta_v_mid'] <- denv3_2014_esimates$beta_v_mid ## seasonality parameters estimated during DENV3-2014 fitting
thetaAlltab[1,iiH,'beta_grad'] <- denv3_2014_esimates$beta_grad
thetaAlltab[1,iiH,'beta_mid'] <- denv3_2014_esimates$beta_mid
thetaAlltab[1,iiH,'beta_width'] <- denv3_2014_esimates$beta_width

if(include.2014.control == T){
  thetaAlltab[1,iiH,'beta_base'] <- denv3_2014_esimates$beta_base
}else{
  thetaAlltab[1,iiH,'beta_base'] <- 0
}

## Plot "control function" with initial starting values
control_plot_vals <- as.data.frame(t(thetaAlltab[1,iiH,]))
control_plot <- sapply(1:max(time.vals), function(xx){control_f(xx, 
                                         base = control_plot_vals$beta_base,
                                         mid = control_plot_vals$beta_mid,
                                         width = control_plot_vals$beta_width
                                         )})
#plot(date.vals[30:100], control_plot[30:100], type="l", ylim=c(0,1), bty="n", ylab="Relative transmission", xlab="Date", xaxt="n")
plot(startdate+(1:max(time.vals))[340:(340+180)], control_plot[340:(340+180)], type="l", ylim=c(0,1), bty="n", ylab="Relative transmission", xlab="Date", xaxt="n")
  axis.Date(side=1, at=seq.Date(startdate+(1:max(time.vals))[340], startdate+(1:max(time.vals))[340+180], by = "1 months"), "months", format = "%b%y")

dev.copy(pdf, here::here("output/fig_supp_controlFn.pdf"), 6,4)
  dev.off()

# run a single simulation with initial values -----------------------------
source(here::here("Rscripts/zika-single-simulation.R"))
initial_beta_h <- thetaAlltab[1,iiH,][["beta_h"]]
thetaAlltab[1,iiH,][["intro_mid"]] <- 650
run1singlesim <- zika_single_sim(initial_beta_h)
if(is.nan(run1singlesim) | run1singlesim==-Inf){
  t_rate <- seq(0.2, 0.5, 0.02)
    lik_search <- sapply(t_rate, function(xx){zika_single_sim(xx)})
    max_beta_h <- t_rate[which(lik_search==max(lik_search, na.rm = T))]
    thetaAlltab[1,iiH,"beta_h"] <- max_beta_h[1]
  zika_single_sim(max_beta_h[1])
}
if(zika_single_sim(thetaAlltab[1,iiH,][["beta_h"]])==-Inf){
  stop("Starting with infinite likelihood - this won't end well!")
}
zika_single_sim(thetaAlltab[1,iiH,][["beta_h"]])

## Plot "introduction function" with initial starting values
intro_plot_vals <- as.data.frame(t(thetaAlltab[1,iiH,]))
intro_plot <- sapply(time.vals, function(xx){intro_f(xx, 
                                                      mid=intro_plot_vals$intro_mid,
                                                      width=intro_plot_vals$intro_width,
                                                      base=intro_plot_vals$intro_base
)})
plot(date.vals[1:30], intro_plot[1:30], type="l", bty="n", ylab="ZIKV introductions", xlab="Date", xaxt="n")
  axis.Date(side=1, at=seq.Date(date.vals[1], date.vals[30], by = "1 months"), "months", format = "%b%y")
integrate(intro_f, -Inf, Inf)

dev.copy(pdf, here::here("output/fig_supp_introFn.pdf"), 6,4)
  dev.off()
  
# Print starting conditions for model run --------------------------
print("Theta initial starting conditions")
c(thetatab[1,], thetaAlltab[1,iiH,])
diag(cov_matrix_thetaAll)

# Run MCMC loop and save results ------------------------------------------
(st.time <- Sys.time())

m = 1; iiM=1
foreach(iiM=multichain) %dopar% {  # Loop over regions with parallel MCMC chains

adapt_size_start <- ceiling(0.01*MCMC.runs)

for (m in 1:MCMC.runs){
  # Scale COV matrices for resampling using error term epsilon
  if(m==1){
    epsilon0 = 0.001
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
    sim_liktab_current = sim_liktab[m]
    prior_current = prior[m]
    covmat.empirical <- cov_matrix_thetaA
    theta.mean <- thetaAlltab_current[1,]
    # initialise counter for storing results (m/thining parameter)
    j=1
  }else{
    #scaling.multiplier <- exp((1-1e-7)^(m-adapt_size_start)*(accept_rate-0.234))
    #epsilon0 <- epsilon0 * scaling.multiplier
    #epsilon0 <- min(epsilon0, 0.5)
    #scaling.multiplier <- exp((0.9999)^(m-adapt_size_start)*(accept_rate-0.234))
    epsilon0 <- max(min(0.1,exp(log(epsilon0)+(accept_rate-0.234)*0.999^m)),1e-6) # Stop epsilon getting too big or small
    
    cov_matrix_theta=epsilon0*cov_matrix_theta0
    cov_matrix_thetaA=epsilon0*cov_matrix_thetaAll
    cov_matrix_theta_init=epsilon0*cov_matrix_theta_initAll
  }
  
  ## DELETED: Resample global theta every 2nd step 
  theta_star <- thetatab_current
  
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
  for(kk in 1){ #itertab){ 
    iiH <- kk
      
    if(m==1){ # Don't resample on 1st step - check the zeroes!
      output_H <- SampleTheta(global=0,thetaAlltab_current[iiH,],theta_initAlltab_current[iiH,],0*cov_matrix_thetaA,0*cov_matrix_theta_init)
    }else{
      output_H <- SampleTheta(global=0,thetaAlltab_current[iiH,],theta_initAlltab_current[iiH,],cov_matrix_thetaA,cov_matrix_theta_init)
    } 
    thetaA_star <- output_H$thetaS
    theta_init_star <- output_H$theta_initS
    
    # Run model simulation
    output1 <- Deterministic_modelR_final_DENVimmmunity_monthly(theta = c(theta_star,thetaA_star,theta_denv), theta_init_star, seroposdates=seroposdates)
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
    prior.star <- (prior.theta$prior.star * prior.star)
    prior.current <- (prior.theta$prior * prior.current)
  } # end loop over regions
  
  # Calculate probability function - MH algorithm
    estimated_params <- colnames(cov_matrix_thetaA)[diag(cov_matrix_thetaA)!=0]
    q_theta_given_theta_star = sum(log(thetaAllstar[iiH, estimated_params]))
    q_theta_star_given_theta = sum(log(thetaAlltab_current[1, estimated_params]))
  
  val <- (prior.star/prior.current) * exp((sim_marg_lik_star - sim_liktab_current) +  (q_theta_given_theta_star - q_theta_star_given_theta)) 
  
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
      sim_liktab[j+1] = sim_liktab[j]
      accepttab[j]=0
      prior[j+1] = prior[j] 
    }
    accept_rate = sum(accepttab[1:j])/j
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
  
  if(m < adapt_size_start){
    accept_rate <- 0.234
  }
  
  if(m %% min(MCMC.runs, 1000)==0){
    print(c(iiM,"m"=m,  
            "accept"=signif(sum(accepttab[1:m])/m,3), 
            "lik"=signif(sim_liktab_current,3),
            "val" = val,
            "prior" = prior.star/prior.current,
            epsilon0
            ),
            digits = 2)
    }  
  if(m %% min(MCMC.runs, 5000)==0){
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
            file=paste("posterior_output/outputR_",iiM,"_",run.name,".RData",sep=""))
  }
  } # End MCMC loop
}
endtime <- Sys.time()
(time.take = endtime-st.time)

# load(here::here("posterior_output/outputR_1_test_monthly.RData"))
# m <-  max(which(!is.na(thetaAlltab[, 1, "rep"])))
# par(mfrow = c(2,2))
# plot(thetaAlltab[floor(m*0.4):m, 1, "beta_h"], type = "l")
# plot(thetaAlltab[floor(m*0.4):m, 1, "rep"], type = "l")
# plot(thetaAlltab[floor(m*0.4):m, 1, "intro_mid"], type = "l")
# plot(thetaAlltab[floor(m*0.4):m, 1, "intro_base"], type = "l")
# plot(c_trace_tab[2,1,], type = 'l')
#   lines(c_trace_tab[m,1,])
#   lines(c_trace_tab[round(m/2),1,])
# plot(accepttab[1:m], type = "p", main = mean(accepttab[1:m], na.rm = T))
# plot(sim_liktab[1:m], type = "l")
# plot(prior[2:m], type = "l")
