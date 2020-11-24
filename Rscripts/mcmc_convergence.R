# - - - - - - - - - - - - - - - - - - - - - - - 
# MCMC convergennce diagnostics
# Author: Alasdair Henderson
# github.com/a-henderson91/fiji-zikv-model
# - - - - - - - - - - - - - - - - - - - - - - - 
pacman::p_load('devtools')
#install_github('sbfnk/fitR')
library('fitR')
pacman::p_load('coda')
pacman::p_load("basicMCMCplots")

virus <- "ZIKV"  # DEN3 or ZIKV

ll <- 1

here::here()
if(virus=="ZIKV"){
  iiH <- 1
  preamblecode <- "preamble_zikvfiji.R"
  post_file_path <- "posterior_output"
}else if(virus=="DEN3"){
  iiH <- 2
  preamblecode <- "preamble_denv3.R"
  post_file_path <- "posterior_denv2014fit"
}
source(here::here("Rscripts", preamblecode))
if(virus == "DEN3"){model1_name <- run.name}

mcmc.burn <- 0.4

# load vectorbornefit code ------------------------------------------------
all_files <- list.files(here::here("Rfunctions"))
sapply(all_files, function(xx){source(here::here("Rfunctions/", xx))})

# load other necessary files ----------------------------------------------
## 
##source(here::here("Rscripts/load_tmrca_calc_prior.R"))

## 
labelN <- 1
iiH <- 2
load_posterior_1 <- load.posteriors(load.run.name=model1_name, file.path="posterior_denv2014fit", iiH, mcmc.burn=mcmc.burn)
  list2env(load_posterior_1,globalenv())
priorFit <- fitdistr(thetatab$beta_base, densfun="normal")
prior_vectorcontrol <- function(x){dnorm(x, mean=priorFit$estimate[1], sd=priorFit$estimate[2])}
iiH <- 1

##
source(here::here("Rscripts/prior_distributions.R"))

# import posteriors -------------------------------------------------------
labelN <- 1
m.tot <- length(list.files(path = here::here(post_file_path), pattern = paste0("*", run.name)))

# load data --------------------------------------------------------------
data <- load.data.multistart(Virus = virusTab[iiH], startdate = start.output.date, serology.file.name = serology.excel, init.values.file.name = init.conditions.excel, add.nulls = 0) #virusTab[iiH], dataTab[iiH])
  list2env(data,globalenv())

# trace plots -------------------------------------------------------------
m.tot <- length(list.files(path = here::here(post_file_path), pattern=paste0("*",run.name)))
# iiM=1
load_theta <- function(iiM, load.run.name, burnin = F, mcmc.burn = NA, load_iiH = 1){
  load(paste0(here::here(post_file_path),"/outputR_",iiM,"_",load.run.name,".RData",sep=""))
  theta_select <- data.frame(thetaAlltab[,load_iiH,])
  
  # Return a character vector of variable names which have 0 variance
  theta_select_names <- names(theta_select)[vapply(theta_select, function(x) var(x, na.rm = T) > 1e-10, logical(1))]
  theta_select_names <- theta_select_names[!is.na(theta_select_names)]
  
  if(burnin == T){
    mcmc_samples <- length(sim_liktab)
    maxB <- sum(sim_liktab!=-Inf)/mcmc_samples
    minB <- mcmc.burn*maxB
    picks <- c(round(minB*mcmc_samples):round(maxB*mcmc_samples))
  }else{
    mcmc_samples <- length(sim_liktab)
    maxB <- sum(sim_liktab!=-Inf)/mcmc_samples
    picks <- c(1:round(maxB*mcmc_samples))
  }
  theta_select <- theta_select[, theta_select_names]
  #mcmc_samples <- length(sim_liktab)
  maxB <- sum(sim_liktab!=-Inf)/mcmc_samples
  
  list(trace = theta_select[picks,], acceptance.rate = sum(accepttab, na.rm = T)/length(accepttab), picks = picks)
}

if(virus == "DEN3"){load_iiH <- 2}else{load_iiH <- 1}
theta_list <- list()
for(iiM in 1:m.tot){
  theta_tmp <- load_theta(iiM, load.run.name = run.name, burnin = F, load_iiH = load_iiH)
  name <- paste0("theta_", iiM)
  theta_list[[name]] <- theta_tmp$trace
}
head(theta_list)

vars_to_plot <- names(theta_list$theta_1)
nn_to_plot <- length(vars_to_plot)
par(mfrow=c(ceiling(nn_to_plot/2),2))

plot_function <- function(var_name, plot_abline = F){
  cols_define <- c(col1, col2, col7)
  plot(theta_list$theta_1[[var_name]], type='l', ylab = labs_to_plot[[var_name]], xlab = "Iteration", col = 0)
  for(iiM in 1:m.tot){
    tmp <- theta_list[[iiM]]
    lines(tmp[[var_name]], type='l', ylab = var_name, xlab = "Iteration", col = cols_define[iiM])
    if(plot_abline == T){
      abline(v = 0.25*length(tmp[[var_name]]), col = 1, lty = 2)
      abline(v = 0.5*length(tmp[[var_name]]), col = 1, lty = 2)
    }
  }
}

# png(here("output/traceplot.png"), 
#      width = 4, height = 6, units = "in",
#      res = 150)
#   par(mfrow=c(ceiling(nn_to_plot/2),2), mar = c(4, 4, 0.5, 0.2))
#   sapply(vars_to_plot, function(xx){plot_function(xx, plot_abline = T)})
# dev.off()

# Define burnin -----------------------------------------------------------
mcmc.burn <- 0.4

# Re-import posteriors with burnin ----------------------------------------
theta_list <- list()
for(iiM in 1:m.tot){
  theta_tmp <- load_theta(iiM, load.run.name = run.name, burnin = T, mcmc.burn = mcmc.burn, load_iiH = load_iiH)
  name <- paste0("theta_", iiM)
  theta_list[[name]] <- theta_tmp$trace
}

if(virus == "ZIKV"){
  labs_to_plot <- c("Transmisssion rate", "Rep. prop.", "Control effect", "Cross-protection prop.",
                    "False Pos.", "Sensitivity",  "intro (mid)", "intro (height)")
  names(labs_to_plot) <- vars_to_plot
}else if(virus == "DEN3"){
  labs_to_plot <- c("Transmisssion rate", "Rep. prop.","Initial immune", "Control effect", "False pos.", "Sensitivity")
  names(labs_to_plot) <- vars_to_plot
}

png(here("output", paste0("traceplot_burnin_", run.name, ".png")), 
     width = 4, height = 6, units = "in",
     res = 150)
  par(mfrow=c(ceiling(nn_to_plot/2),2), mar = c(4, 4, 0.5, 0.2))
  sapply(vars_to_plot, function(xx){plot_function(xx)})
dev.off()

# histograms of output ----------------------------------------------------
## combine the chains
theta_joint <- rbind(theta_list$theta_1, theta_list$theta_2)

if(virus == "DEN3"){
prior_vectorcontrol <- function(x){
  dnorm(x, mean = 0.57, sd = 0.15) ## based on Kucharski et al. 2018 Elife
}
priorRec0 <- function(x){dnorm(x, 0.331, 0.25)}

priorIntro <- function(x){dunif(x, min=0, max=307)} ## DENV import mid point is 2013-11-07 (first case reported) at latest
priorInitInf <- function(x){dunif(x, min=0, max=1000)} ## DENV import max is 1000 in one day
priorEpsilon <- function(x){dtruncnorm(x, a=0, b=1, mean=0, sd=0.15)} ## 100% sensitivity and specificity in control panel 
priorAlpha <- function(x){dtruncnorm(x, a=0, b=1, mean=1, sd=0.15)}
}

pdf(here("output", paste0("mcmc_density_plots_", run.name, ".pdf")), height = 6, width = 6)
par(mfrow=c(ceiling(nn_to_plot/2),2), mar = c(4, 4, 2, 4))
## beta, beta_v, beta_base, rep, chi, epsilon, rho, intro_mid, intro_base, intro_width
hist(theta_joint$beta_h, main = "", xlab = "") 
par(new=T)
curve(priorBeta, col = 2, lwd = 2, yaxt = "n", xaxt = "n", main = "", ylab = "", xlab = "")
axis(side = 4)
mtext(side = 4, "Prior", cex = 0.7, padj = 3)
mtext("Transmission rate", side=3, adj=0, font=2)

hist(theta_joint$rep, main = "", xlab = "") 
mtext("Reporting prop.", side=3, adj=0, font=2)
  
hist(theta_joint$beta_base, main = "", xlim = c(0,1), xlab = "") 
par(new=T)
curve(prior_vectorcontrol, col = 2, lwd = 2, yaxt = "n", xaxt = "n", main = "", ylab = "", xlab = "", from = 0, to = 1)
axis(side = 4)
mtext(side = 4, "Prior", cex = 0.7, padj = 3)
mtext("Control effect", side=3, adj=0, font=2)

if(virus == "DEN3"){
  hist(theta_joint$rec0, main = "", xlim = c(0,1), xlab = "") 
  par(new=T)
  curve(priorRec0, col = 2, lwd = 2, yaxt = "n", xaxt = "n", main = "", ylab = "", xlab = "", from = 0, to = 1)
  axis(side = 4)
  mtext(side = 4, "Prior", cex = 0.7, padj = 3)
  mtext("Initial immune", side=3, adj=0, font=2)
  
  
}

if(virus == "ZIKV"){
hist(theta_joint$epsilon, main = "", xlab = "", xlim = c(0,1)) 
par(new=T)
curve(priorEpsilon, col = 2, lwd = 2, yaxt = "n", xaxt = "n", main = "", ylab = "", xlab = "", xlim = c(0,1))
axis(side = 4)
mtext(side = 4, "Prior", cex = 0.7, padj = 3)
mtext("False positivity", side=3, adj=0, font=2)

hist(theta_joint$alpha, main = "", xlab = "", xlim = c(0,1)) 
par(new=T)
curve(priorAlpha, col = 2, lwd = 2, yaxt = "n", xaxt = "n", main = "", ylab = "", xlab = "")
axis(side = 4)
mtext(side = 4, "Prior", cex = 0.7, padj = 3)
mtext("Sensitivity", side=3, adj=0, font=2)

  
  hist(theta_joint$chi, main = "", xlab = "") 
  par(new=T)
  curve(priorChi(x), col = 2, lwd = 2, yaxt = "n", xaxt = "n", main = "", ylab = "", xlab = "")
  axis(side = 4)
  mtext(side = 4, "Prior", cex = 0.7, padj = 3)
  mtext("Cross-protect", side=3, adj=0, font=2)

  hist(theta_joint$intro_mid, main = "", xlim = c(0,1500), xlab = "")
  par(new=T)
  curve(priorIntro, col = 2, lwd = 2, yaxt = "n", xaxt = "n", main = "", ylab = "", xlab = "", from = 0, to = 1500)
  axis(side = 4)
  mtext(side = 4, "Prior", cex = 0.7, padj = 3)
  mtext("Intro mid", side=3, adj=0, font=2)
  
  hist(theta_joint$intro_base, main = "", xlim = c(0,10), xlab = "") 
  par(new=T)
  curve(priorInitInf, col = 2, lwd = 2, yaxt = "n", xaxt = "n", main = "", ylab = "", xlab = "", from = 0, to = 10)
  axis(side = 4)
  mtext(side = 4, "Prior", cex = 0.7, padj = 3)
  mtext("Intro height", side=3, adj=0, font=2)
}

dev.off()

# correlation plots -------------------------------------------------------
sample.p <- sample(length(theta_joint$beta_h),100,replace=T)
thinner.theta <- theta_joint[sample.p,]

png(here("output", paste0("corrplot_", run.name, ".png")), 
     width = 12, height = 12, units = "in",
     res = 200)
par(mfrow = c(length(vars_to_plot), length(vars_to_plot)),  mar = c(3, 4, 2, 1) )
for(ii in 1:length(vars_to_plot)){
  for(jj in 1:length(vars_to_plot)){
    if(ii<=jj){
      if(ii == jj){
        hist(theta_joint[[vars_to_plot[ii]]],xlab=labs_to_plot[ii],main=NULL)
        mtext(labs_to_plot[ii], adj = 1, font = 1, cex = 0.9)
      }else{
        plot(thinner.theta[[vars_to_plot[ii]]],thinner.theta[[vars_to_plot[jj]]],pch=19,cex=0.2, xlab="", ylab="",col="white",xaxt="n",yaxt="n",axes=F)
      }
    }else{
      plot(thinner.theta[[vars_to_plot[ii]]],thinner.theta[[vars_to_plot[jj]]],pch=19,cex=0.3, xlab=vars_to_plot[ii], ylab=vars_to_plot[jj])
      points(median(thinner.theta[[vars_to_plot[ii]]]),median(thinner.theta[[vars_to_plot[jj]]]),col="orange",cex=1.5,lwd=2)
    }
  }
}
dev.off()

# Diagnostics table -------------------------------------------------------
ess_theta <- apply(theta_joint, 2, function(xx){effectiveSize(xx)})
write.csv(signif(ess_theta,3), here::here("output", "ess_theta.csv"))

#x_ess_theta <- xtable::xtable(t(ess_theta), include.rownames = F)  
#xtable::print.xtable(x_ess_theta, include.rownames = F)
