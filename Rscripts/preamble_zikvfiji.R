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
#seasonal.transmission   <- T # whether to estimate seasonal transmission or not. If false - zeroes BETA_V_AMP after estimation in main script
include.sero.likelihood <- T
#vector.control          <- T # 
include.2014.control    <- T # if False then beta_base set to 0
run.name <- "0219_github" 
## MCMC parameters 
MCMC.runs <- 100 #number of MCMC iterations 
thinning.parameter <- 1
multichain <- c(1:3)  # n chains to run in parallel
mcmc.burn <- 0.4

# data file names ---------------------------------------------------------
locationtab <- c("Zika2016","Central2014") 
dataTab <- c("Central_Fiji_2016Z_timeseries","Central_2014_timeseries") 
virusTab <- c("Zika","DEN3")
serology.excel <- "Fiji_serology"
init.conditions.excel <- "thetaR_IC_zika"
iiH <- 1
locnn <- 2 ## needs to be >1 to set up results matrix not vector
itertab <- c(1:locnn)
itertabM=c(1:1) # Iterate over locations in set up and MCMC


# model parameters --------------------------------------------------------
startdate <- as.Date("2013-01-01")  ## shouldn't be changed (some parameters are relative to this date)
dt <-  7                            # time step
start.output.date <- as.Date("2013-01-01")
end.output.date <- as.Date("2018-01-01")
seroposdates <- c(as.Date("2013-11-15"), as.Date("2015-11-02"), as.Date("2017-06-18")) ## dates of seroprevalence surveys

# print model run info ------------------------------------------------
print(paste0("sample start time? ", sample.start.point))
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
