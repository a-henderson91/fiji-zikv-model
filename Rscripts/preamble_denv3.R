# - - - - - - - - - - - - - - - - - - - - - - - 
# Transmission modelling preamble
# Author: Alasdair Henderson
# github.com/a-henderson91/fiji-zikv-model
# - - - - - - - - - - - - - - - - - - - - - - - 

# Load packages and set up pathogen -----------------------------------------------
library(coda)
library(colorspace)
library(doMC)
clust1<-registerDoMC(4)  #change to your number of CPU cores
library(data.table)
library(deSolve)
library(foreach)
library(ggplot2)
library(gridExtra)
library(lubridate)
library(magrittr)
library(MASS)
library(mvtnorm)
library(tidyverse)
library(truncnorm)

library(here)

# model run options to change ---------------------------------------------
seasonal.transmission   <- T # whether to estimate seasonal transmission or not. If false - zeroes BETA_V_AMP after estimation in main script
include.sero.likelihood <- T
vector.control          <- T # 
include.2014.control    <- T # if False then beta_base set to 0
run.name <- "0219_model1_denv3" 

## MCMC parameters 
MCMC.runs <- 1e4 #number of MCMC iterations 
thinning.parameter <- 1
multichain <- c(1:3)  # n chains to run in parallel
mcmc.burn <- 0.4

# data file names ---------------------------------------------------------
#locationtab <- c("Central2014","Zika2016") # IF running Central2014 - this needs to go in position no.1
#dataTab <- c("Central_2014_timeseries","Central_Fiji_2016Z_timeseries") 
#virusTab <- c("DEN3","Zika")
locationtab <- c("Zika2016","Central2014") # IF running Central2014 - this needs to go in position no.1
dataTab <- c("Central_Fiji_2016Z_timeseries","Central_2014_timeseries") 
virusTab <- c("Zika","DEN3")
serology.excel <- "Fiji_serology"
init.conditions.excel <- "thetaR_IC_denv3"
#iiH <- 1
locnn <- 2 ## needs to be >1 to set up results matrix not vector
itertab <- c(1:locnn)
itertabM=c(1:1) # Iterate over locations in set up and MCMC

# model parameters --------------------------------------------------------
startdate <- as.Date("2013-01-01")  ## shouldn't be changed (some parameters are relative to this date)
dt <-  7                            # time step
start.output.date <- as.Date("2013-01-01")
end.output.date <- as.Date("2015-11-02")
seroposdates <- c(as.Date("2013-11-15"), as.Date("2015-11-02"))

# print model run info ------------------------------------------------
print(paste0("MCMC runs = ", MCMC.runs))
print(paste0("start date = ", as.Date(startdate)))
print(paste0("No. serology dates = ", length(seroposdates)))
print(paste0("Run name = ", run.name))
print(paste0("No. of chains = ", length(multichain)))

# colours for plots -------------------------------------------------------
datacol=rgb(0.4,0.4,0.4)
col2=rgb(0.8,0,0.5,0.8)
col2a=rgb(0.8,0,0.5,0.4)
col1=rgb(0,0.3,1,0.8)
col1a=rgb(0,0.3,1,0.4)

col3=rgb(0.4,0.4,0.4,0.8)
col3a=rgb(0.4,0.4,0.4,0.4)
colgrid=rgb(0.4,0.4,0.4,0.1)
col4=rgb(1,0.7,0.2,0.8)
col4a=rgb(1,0.7,0.2,0.4)

col5=rgb(1,0.4,0.1,0.8)
col5a=rgb(1,0.4,0.1,0.4)
col6=rgb(0.9,0.1,0.1,0.8)
col6a=rgb(0.9,0.1,0.1,0.4)

col7=rgb(0.02,0.3,0.02,0.8)
col7a=rgb(0.02,0.3,0.02,0.4)


col1 <- "#3880bf"
col2 <- "#ef286e"
col3 <- "#c06b74"
col4 <- "#3eeaef"
col5 <- "#881448"
col6 <- "#0ba47e"
col7 <- "#1e5c4a"

col1a <- alpha("#3880bf" ,0.4)
col2a <- alpha("#ef286e" ,0.4)
col3a <- alpha("#c06b74" ,0.4)
col4a <- alpha("#3eeaef" ,0.4)
col5a <- alpha("#881448" ,0.4)
col6a <- alpha("#0ba47e" ,0.4)
col7a <- alpha("#1e5c4a" ,0.4)

theme_zika_fiji <- function(){
  theme(title = element_text(size = 12, face = "bold"),
        panel.background = element_rect(fill = "white",
                                        colour = "white",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(size=.1, color=rgb(0,0,0,0.3) ) ,
        axis.ticks = element_blank(),
        line = element_line(colour = "gray70"),
        rect = element_rect(colour = "black", fill = "white")
  )
}









# Fit FP Zika epidemic ----------------------------------------------------
#zikaFP <- read.csv("data/FP_zika.csv", stringsAsFactors = F)
#dd <- substr(zikaFP$date,1,2)
#mm <- substr(zikaFP$date,4,5)
#Y <- substr(zikaFP$date,7,10)
#zikaFP$realdate <- paste0(Y,"-",mm,"-",dd)
#zikaFP$realdate <- as.Date(zikaFP$realdate)
#
#zikaFPfull <- data.frame("date"=zikaFP$realdate)
#full_zika <- with(zikaFP,zikaFP$tahiti + zikaFP$ile_sous + zikaFP$moorea + zikaFP$tuamotu + zikaFP$marquises + zikaFP$australes)
#avg_rep <- mean(c(zikaFP$tahitiRep,zikaFP$ile_sousRep,zikaFP$mooreaRep,zikaFP$tuamotuRep,zikaFP$marquisesRep,zikaFP$australesRep))
#
#full_zika_raw <- full_zika
#full_zika <- full_zika*(1/avg_rep)
#
#zikaFPfull <- cbind(zikaFPfull, "rawZIKV"=full_zika_raw, "ZIKV"=full_zika)
#zikaFPfull$time <- seq(1,length(zikaFPfull$date)*7,7)
#
##zikaFPfull <- zikaFPfull[1:18,]
#quadraticModel <- lm(zikaFPfull$ZIKV ~ zikaFPfull$time + I(zikaFPfull$time^2))
#x<-sort(zikaFPfull$time)
#y<-quadraticModel$fitted.values[order(zikaFPfull$time)]
#m <- quadraticModel$coefficients[1]
#beta1 <- quadraticModel$coefficients[2]
#beta2 <- quadraticModel$coefficients[3]
#x <- zikaFPfull$time
#y <- m + beta1*x + beta2*(x^2)
#
#offset=0#-(as.numeric(NEW.DATE-as.Date("2013-10-11")))
#
#Ctreg <- function(t, offset){
#  if(offset+t<0){
#    y2 <- 0
#    as.numeric(y2)
#  }else{
#    tt <- offset+t
#    y <- m + beta1*tt + beta2*(tt^2)
#    y2 <- max(y, 0)
#    as.numeric(y2)#*10) #reporting rate = 10%
#  }
#}

# Fit DENV3 epidemic ------------------------------------------------------
#denv3 <- read.csv("data/Central_2014_timeseries.csv", stringsAsFactors = F)
#
#denv3$time <- seq(1,length(denv3$date)*7,7)
#
#denv3 <- denv3[,c(1:2,6)]
#colnames(denv3)[2] <- "DENV3"
#poissonModel <- glm(denv3$DENV3 ~ denv3$time + I(denv3$time^2), family="poisson")
#
#x<-sort(denv3$time)
#y<-poissonModel$fitted.values[order(denv3$time)]
#d.m <- poissonModel$coefficients[1]
#d.beta1 <- poissonModel$coefficients[2]
#d.beta2 <- poissonModel$coefficients[3]
#d.x <- denv3$time
#d.y = exp(d.m + d.beta1*d.x + d.beta2*(d.x^2))
#
#D3reg <- function(t){
#  #y <- exp(d.m + d.beta1*t + d.beta2*(t^2))
#  y <- exp(1.58971 + 0.1019994*t - 0.0005130412*(t^2))
#  y2 <- max(y, 0)
#  as.numeric(y2*(1/0.065)) #reporting rate = 10%
#}
#
#DENVinfections <- as.numeric(y*(1/0.1))
