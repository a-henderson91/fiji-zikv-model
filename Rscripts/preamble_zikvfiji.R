# - - - - - - - - - - - - - - - - - - - - - - - 
# Transmission modelling preamble for ZIKV
# Author: Alasdair Henderson
# github.com/a-henderson91/fiji-zikv-model
# - - - - - - - - - - - - - - - - - - - - - - - 

# Load packages and set up pathogen -----------------------------------------------
#install.packages("pacman")
pacman::p_load(coda,
               here,
               lubridate,
               colorspace,
               doMC,
               data.table,
               deSolve,
               foreach,
               ggplot2,
               gridExtra,
               magrittr,
               MASS,
               mvtnorm,
               stringr,
               tidyverse,
               truncnorm)
library(here)
clust1 <- registerDoMC(4)  #change to your number of CPU cores

# model run options to change ---------------------------------------------
seasonal.transmission   <- T # whether to estimate seasonal tra nsmission or not. If false - zeroes BETA_V_AMP after estimation in main script
include.sero.likelihood <- T # whether to include serological data in likelihood
vector.control          <- T # whether to reduce transmission in March2014 when vector control campaign was in effect 
include.2014.control    <- T # if False then beta_base set to 0
limit.to.2013           <- F # if True then prior on intro time is limited to 365 (i.e. ZIKV starts in 2013)

run.name <- "1123_mainZIKV"
model1_name <- "1123_model1"
dt <- (7*52)/12
  
## MCMC parameters 
MCMC.runs <- 2e5 #number of MCMC iterations 
thinning.parameter <- 1
multichain <- c(1:3)  # n chains to run in parallel
mcmc.burn <- 0.4

# data file names ---------------------------------------------------------
locationtab <- c("Zika2016","Central2014") 
dataTab <- c("Central_Fiji_2016Z_timeseries","Central_2014_timeseries") 
virusTab <- c("ZIKV","DEN3")
serology.excel <- "Fiji_serology"
init.conditions.excel <- "thetaR_IC_zika"
#iiH <- 1
locnn <- 2 ## needs to be >1 to set up results matrix not vector
itertab <- c(1:locnn)
itertabM=c(1:1) # Iterate over locations in set up and MCMC

# model parameters --------------------------------------------------------
startdate <- as.Date("2013-01-01")  ## shouldn't be changed (some parameters are relative to this date)
start.output.date <- as.Date("2013-01-01")
end.output.date <- as.Date("2018-01-01")
seroposdates <- c(as.Date("2013-11-15"), as.Date("2015-11-02"), as.Date("2017-06-18")) ## dates of seroprevalence surveys

# print model run info ------------------------------------------------
print(paste0("MCMC runs = ", MCMC.runs))
print(paste0("start date = ", as.Date(startdate)))
print(paste0("No. serology dates = ", length(seroposdates)))
print(paste0("Run name = ", run.name))
print(paste0("No. of chains = ", length(multichain)))
print(paste0("Name of model1 results = ", model1_name))
print(paste0("Name of run = ", run.name))

# colours for plots -------------------------------------------------------
datacol=rgb(0.4,0.4,0.4)
col1=rgb(0,0.3,1,0.8)
col1a=rgb(0,0.3,1,0.4)
col2=rgb(0.8,0,0.5,0.8)
col2a=rgb(0.8,0,0.5,0.4)

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


#col1 <- "#3880bf"
#col2 <- "#ef286e"
#col3 <- "#c06b74"
#col4 <- "#3eeaef"
#col5 <- "#881448"
#col6 <- "#0ba47e"
#col7 <- "#1e5c4a"
#
#col1a <- alpha("#3880bf" ,0.4)
#col2a <- alpha("#ef286e" ,0.4)
#col3a <- alpha("#c06b74" ,0.4)
#col4a <- alpha("#3eeaef" ,0.4)
#col5a <- alpha("#881448" ,0.4)
#col6a <- alpha("#0ba47e" ,0.4)
#col7a <- alpha("#1e5c4a" ,0.4)

theme_zika_fiji <- function(){
  theme(plot.title = element_text(size = 12, hjust=0),
        axis.title = element_text(colour = rgb(0,0,0,0.8)), 
        panel.background = element_rect(fill = "white",
                                        colour = "white",
                                        size = 0.5, linetype = "solid"),
        axis.line = element_line(size=rel(0.5), colour = rgb(0,0,0,0.8)),
        panel.grid.major = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_line(size=.2, color=rgb(0,0,0,0.2)) ,
        panel.grid.major.y = element_line(size=.2, color=rgb(0,0,0,0.3)) ,
        axis.ticks = element_line(size=.1, color=rgb(0,0,0,0.4)),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = rel(0.75), hjust = 0),
        ##line = element_line(colour = "gray70"),
        ##rect = element_rect(fill = "white"),
        legend.background = element_blank(), legend.key = element_blank()
  )
}
