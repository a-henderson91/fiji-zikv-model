# - - - - - - - - - - - - - - - - - - - - - - - 
# MCMC convergennce diagnostics
# Author: Alasdair Henderson
# github.com/a-henderson91/fiji-zikv-model
# - - - - - - - - - - - - - - - - - - - - - - - 
library('devtools')
#install_github('sbfnk/fitR')
library('fitR')
library('coda')
library("basicMCMCplots")

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

mcmc.burn <- 0.5

# load vectorbornefit code ------------------------------------------------
all_files <- list.files(here::here("Rfunctions"))
sapply(all_files, function(xx){source(here::here("Rfunctions/", xx))})

# load other necessary files ----------------------------------------------
## 
source(here::here("Rscripts/load_tmrca_calc_prior.R"))

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
load_theta <- function(iiM, load.run.name, burnin = F, mcmc.burn = NA){
  load(paste0(here::here(post_file_path),"/outputR_",iiM,"_",load.run.name,".RData",sep=""))
  theta_select <- data.frame(thetaAlltab[,iiH,])
  
  # Return a character vector of variable names which have 0 variance
  theta_select_names <- names(theta_select)[vapply(theta_select, function(x) var(x, na.rm = T) > 1e-10, logical(1))]
  
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

theta_list <- list()
for(iiM in 1:m.tot){
  theta_tmp <- load_theta(iiM, load.run.name = run.name, burnin = F)
  name <- paste0("theta_", iiM)
  theta_list[[name]] <- theta_tmp$trace
}
head(theta_list)

vars_to_plot <- names(theta_list$theta_1)
nn_to_plot <- length(vars_to_plot)
par(mfrow=c(ceiling(nn_to_plot/2),2))

plot_function <- function(var_name, plot_abline = F){
  cols_define <- c(col1, col2)
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

png(here("output/traceplot.png"), 
     width = 4, height = 6, units = "in",
     res = 150)
  par(mfrow=c(ceiling(nn_to_plot/2),2), mar = c(4, 4, 0.5, 0.2))
  sapply(vars_to_plot, function(xx){plot_function(xx, plot_abline = T)})
dev.off()

# Define burnin -----------------------------------------------------------
mcmc.burn <- 0.5

# Re-import posteriors with burnin ----------------------------------------
theta_list <- list()
for(iiM in 1:m.tot){
  theta_tmp <- load_theta(iiM, load.run.name = run.name, burnin = T, mcmc.burn = mcmc.burn)
  name <- paste0("theta_", iiM)
  theta_list[[name]] <- theta_tmp$trace
}
head(theta_list)

labs_to_plot <- c("beta (humans)", "beta (mosquitoes)", "rep. prop.", "control", "cross-protect",
                  "assay spec.", "AB duration", "intro (mid)", "intro (height)", "intro (width)")
names(labs_to_plot) <- vars_to_plot

png(here("output/traceplot_burnin.png"), 
     width = 4, height = 6, units = "in",
     res = 150)
par(mfrow=c(ceiling(nn_to_plot/2),2), mar = c(4, 4, 0.5, 0.2))
sapply(vars_to_plot, function(xx){plot_function(xx)})
dev.off()

# histograms of output ----------------------------------------------------
## combine the chains
theta_joint <- rbind(theta_list$theta_1, theta_list$theta_2)

pdf(here::here("output", "mcmc_density_plots.pdf"), height = 6, width = 6)
par(mfrow=c(ceiling(nn_to_plot/2),2), mar = c(4, 4, 2, 4))
## beta, beta_v, beta_base, rep, chi, epsilon, rho, intro_mid, intro_base, intro_width
hist(theta_joint$beta_h, main = "", xlab = "") 
par(new=T)
curve(priorBeta, col = 2, lwd = 2, yaxt = "n", xaxt = "n", main = "", ylab = "", xlab = "")
axis(side = 4)
mtext(side = 4, "Prior", cex = 0.7, padj = 3)
mtext("beta (humans)", side=3, adj=0, font=2)

hist(theta_joint$beta_v, main = "", xlab = "") 
par(new=T)
curve(priorBetaV(x), col = 2, lwd = 2, yaxt = "n", xaxt = "n", main = "", ylab = "", xlab = "")
axis(side = 4)
mtext(side = 4, "Prior", cex = 0.7, padj = 3)
mtext("beta (vector)", side=3, adj=0, font=2)

hist(theta_joint$rep, main = "", xlab = "") 
mtext("rep. rate", side=3, adj=0, font=2)

hist(theta_joint$beta_base, main = "") 
par(new=T)
curve(prior_vectorcontrol(x), col = 2, lwd = 2, yaxt = "n", xaxt = "n", main = "", ylab = "", xlab = "")
axis(side = 4)
mtext(side = 4, "Prior", cex = 0.7, padj = 3)
mtext("control", side=3, adj=0, font=2)

hist(theta_joint$chi, main = "", xlab = "") 
par(new=T)
curve(priorChi(x), col = 2, lwd = 2, yaxt = "n", xaxt = "n", main = "", ylab = "", xlab = "")
axis(side = 4)
mtext(side = 4, "Prior", cex = 0.7, padj = 3)
mtext("cross-protect", side=3, adj=0, font=2)

hist(theta_joint$epsilon, main = "", xlab = "") 
par(new=T)
curve(priorepsilon, col = 2, lwd = 2, yaxt = "n", xaxt = "n", main = "", ylab = "", xlab = "")
axis(side = 4)
mtext(side = 4, "Prior", cex = 0.7, padj = 3)
mtext("assay spec.", side=3, adj=0, font=2)

hist(theta_joint$rho, main = "", xlab = "") 
mtext("AB duration", side=3, adj=0, font=2)

hist(theta_joint$intro_mid, main = "", xlim = c(0,1500), xlab = "")
par(new=T)
curve(priorIntro, col = 2, lwd = 2, yaxt = "n", xaxt = "n", main = "", ylab = "", xlab = "", from = 0, to = 1500)
axis(side = 4)
mtext(side = 4, "Prior", cex = 0.7, padj = 3)
mtext("intro mid", side=3, adj=0, font=2)

hist(theta_joint$intro_base, main = "", xlim = c(0,10), xlab = "") 
par(new=T)
curve(priorInitInf, col = 2, lwd = 2, yaxt = "n", xaxt = "n", main = "", ylab = "", xlab = "", from = 0, to = 10)
axis(side = 4)
mtext(side = 4, "Prior", cex = 0.7, padj = 3)
mtext("intro height", side=3, adj=0, font=2)

hist(theta_joint$intro_width, main = "", xlim = c(0,60), xlab = "") 
par(new=T)
curve(priorInitWidth, col = 2, lwd = 2, yaxt = "n", xaxt = "n", main = "", ylab = "", xlab = "", from = 0, to = 60)
axis(side = 4)
mtext(side = 4, "Prior", cex = 0.7, padj = 3)
mtext("intro width", side=3, adj=0, font=2)

dev.off()

# correlation plots -------------------------------------------------------
sample.p <- sample(length(theta_joint$beta_h),1000,replace=T)
thinner.theta <- theta_joint[sample.p,]

png(here("output/corrplot.png"), 
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

## just plot Betas against each other 
par(mfrow = c(2,2),  mar = c(3, 4, 2, 1))
## beta_h hist
hist(theta_joint[[vars_to_plot[1]]],xlab=labs_to_plot[1],main=NULL)
mtext(labs_to_plot[1], adj = 0, font = 1, cex = 0.9)
## blank
plot(thinner.theta[[vars_to_plot[1]]],thinner.theta[[vars_to_plot[1]]],pch=19,cex=0.2, xlab="", ylab="",col="white",xaxt="n",yaxt="n",axes=F)
## betas correlation plot
plot(thinner.theta[[vars_to_plot[1]]],thinner.theta[[vars_to_plot[2]]],pch=19,cex=0.3, xlab=labs_to_plot[1], ylab=labs_to_plot[2])
points(median(thinner.theta[[vars_to_plot[1]]]),median(thinner.theta[[vars_to_plot[2]]]),col="orange",cex=1.5,lwd=2)
## beta_v hist
hist(theta_joint[[vars_to_plot[2]]],xlab=labs_to_plot[2],main=NULL)
mtext(labs_to_plot[2], adj = 0.5, font = 1, cex = 0.9)
dev.copy(pdf, here::here("outputs", "betas_corr.pdf"), 6, 4)
  dev.off()
# Diagnostics table -------------------------------------------------------
ess_theta <- apply(theta_joint, 2, function(xx){effectiveSize(xx)})
write.csv(signif(ess_theta,3), here::here("output", "ess_theta.csv"))

x_ess_theta <- xtable::xtable(t(ess_theta), include.rownames = F)  
xtable::print.xtable(x_ess_theta, include.rownames = F)
