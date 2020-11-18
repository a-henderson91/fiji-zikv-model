##' load.posteriors
##' 
##' Load and combine posterior chains following MCMC run
#' @param agestructure Binary indicator variable if model is age structured (1) or not (0) between children and adults. Defaults to NULL in which case no sampling happens 
#' @param iiH Locationtab index
#' @param mcmc.burn Between 0 and 1. Proportion of iterations to discard.
##' @export 
##' 
#load.run.name=run.name; file.path="posterior_output"; iiH; mcmc.burn=mcmc.burn

load.posteriors <- function(load.run.name, file.path="posterior_outputZ", iiH, mcmc.burn){
thetatabA=NULL
m.tot <- length(list.files(path = paste0(file.path,"/"), pattern = paste0("*",load.run.name)))
theta_inittabA=NULL
c_trace_tab0=NULL
cd_trace_tab0=NULL
s_trace_tab0=NULL
r_trace_tab0=NULL
sim_liktab0=NULL

# iiM=1
for(iiM in 1:m.tot){
  load(paste0(file.path,"/outputR_",iiM,"_",load.run.name,".RData",sep=""))

  thetatab=cbind(data.frame(thetatab),data.frame(thetaAlltab[,iiH,]))
  theta_inittab=data.frame(theta_initAlltab[,iiH,])
  
  mcmc_samples=length(sim_liktab)
  maxB=sum(sim_liktab!=-Inf)/mcmc_samples
  minB=mcmc.burn*maxB
  picks=c(round(minB*mcmc_samples):round(maxB*mcmc_samples))
  
  thetatabA=rbind(thetatabA,thetatab[picks,])
  theta_inittabA=rbind(theta_inittabA,theta_inittab[picks,])
  
  c_trace_tab0 = rbind(c_trace_tab0,c_trace_tab[picks,iiH,])
  if(!is.null(cd_trace_tab0)){
    cd_trace_tab0 = rbind(cd_trace_tab0,cd_trace_tab[picks,iiH,])
  }
  s_trace_tab0 = rbind(s_trace_tab0,s_trace_tab[picks,iiH,])
  r_trace_tab0 = rbind(r_trace_tab0,r_trace_tab[picks,iiH,])
  sim_liktab0  = cbind(sim_liktab0,sim_liktab[picks])
  
}

picks=c(1:length(thetatabA[,1]))
thetatab=thetatabA
theta_inittab=theta_inittabA
c_trace_tab=c_trace_tab0
cd_trace_tab=cd_trace_tab0
s_trace_tab=s_trace_tab0
r_trace_tab=r_trace_tab0
sim_liktab=sim_liktab0

return(list(picks=picks,
         thetatab=thetatab,
         theta_inittab=theta_inittab,
         c_trace_tab=c_trace_tab,
         cd_trace_tab=cd_trace_tab,
         s_trace_tab=s_trace_tab,
         r_trace_tab=r_trace_tab,
         sim_liktab=sim_liktab,
         accepttab=accepttab
         ))
}
