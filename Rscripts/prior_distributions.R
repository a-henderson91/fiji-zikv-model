# - - - - - - - - - - - - - - - - - - - - - - - 
# Prior distibutions for vector borne disease model  
# Author: Alasdair Henderson
# github.com/a-henderson91/fiji-zikv-model
# - - - - - - - - - - - - - - - - - - - - - - - 

var_prior <- 0.1

if(virusTab[1]=="DEN3"){
  prior_p_VEx <- c(10,var_prior)
  prior_p_Exp <- c(5.9,var_prior)
  prior_p_MuV <- c(8,var_prior)
  prior_p_Inf <- c(5,var_prior)
}else{
  prior_p_VEx <- c(15,var_prior)
  prior_p_Exp <- c(5.9,var_prior)
  prior_p_MuV <- c(8.5,var_prior)
  prior_p_Inf <- c(5,var_prior)
}

priorExp <- function(x){dgamma(x,shape=prior_p_Exp[1]/(prior_p_Exp[2]), scale=prior_p_Exp[2])} #1+0*x} #
priorInf <- function(x){dgamma(x,shape=prior_p_Inf[1]/(prior_p_Inf[2]), scale=prior_p_Inf[2])} #1+0*x} #
priorVEx <- function(x){dgamma(x,shape=prior_p_VEx[1]/(prior_p_VEx[2]), scale=prior_p_VEx[2])} #1+0*x} #
priorMuV <- function(x){dgamma(x,shape=prior_p_MuV[1]/(prior_p_MuV[2]), scale=prior_p_MuV[2])} #1+0*x} #
priorOmega <- function(x){dtruncnorm(x,a=1, b=Inf, mean=90, sd=10)} 
priorRec0  <- function(x){dbinom(round(x*nPOP[1]), size=nPOP[1], prob=nLUM[1]/nPOP[1])} 

priorBeta <- function(x){dunif(x, min = 0, max = 1)} 

priorIntro <- function(x){dnorm(x, mean=intro_prior_mu, sd=intro_prior_sigma)} ## from 'Export start time Central Division.R'

priorInitInf <- function(x){dunif(x, min = 1, max = 20)} ## peaks at >=1, no more than 1000 (total into = 4x10xInitInf)
#priorInitWidth <- function(x){dunif(x, min=0, max=50)} ## removed estimate of width

priorChi<-function(x){dunif(x,min=0,max=1)} 
priorepsilon <- function(x){dunif(x, min=0, max=1)}

priorBeta_amp<-function(x){dunif(x,min=0, max=10)} 
priorBeta_mid<-function(x){dunif(x,min=0, max=10)} 