#' Setup results vectors and dataframes and load initial values. Totally susceptible initially
#' 
#' Setup results vectors and dataframes and load initial values. Totally susceptible initially
#' @param iiH PLace in locationtab of location being processed
#' @param parameter_est_file Name of csv file with parameters to be estimated included
#' @export

results_set_up <- function(iiH, parameter_est_file){
  # Set up vectors for storing their values during MCMC loop
  thetaAll=data.frame(rep(NA,locnn))
  for (i in 1:(length(thetaR_IC_local$param)-1)){
    thetaAll <- cbind(thetaAll, rep(NA,locnn))  
  }
  names(thetaAll) <- thetaR_IC_local$param

  compartments=c('s_init','e_init','i1_init','r_init',
                 'sd_init','ed_init','id_init','t1d_init','t2d_init',
                 'sm_init','em_init','im_init',
                 'fpS_init','fpI_init')   
  
  theta_initAll=data.frame(rep(NA,locnn))
  for (i in 1:(length(compartments)-1)){
    theta_initAll <- cbind(theta_initAll, rep(NA,locnn))  
  }
  names(theta_initAll) <- compartments


# Loop to load data for specific location & set initial conditions
for(iiH in itertab){
  ##  Global parameters
  c1=(names(thetaR_IC_local)==locationtab[iiH])
  
  theta <- NULL
  for (i in thetaR_IC_global$param){
    theta <- c(theta, thetaR_IC_global[thetaR_IC_global$param==i, c1])  
  }
  names(theta) <- thetaR_IC_global$param
  
  # Scale age groups by total size
  popsize=theta["npop"]
  popsizeC=theta["npopC"]
  popsizeA=theta["npopA"]
  popsizeTot=popsize #(popsizeC+popsizeA); names(popsizeTot)='npop'
  popsizeC=round(popsize* popsizeC/popsizeTot); popsizeA=round(popsize* popsizeA/popsizeTot)
  
  ## Local parameters
  for (i in thetaR_IC_local$param){
    thetaAll[iiH, names(thetaAll)==i] <- thetaR_IC_local[thetaR_IC_local$param==i, c1] 
  }
  
  if(sum(names(thetaAll)=="Exp")>0){
    thetaAll[iiH,"Exp"] <- 1/thetaAll[iiH,"Exp"]
  }
  if(sum(names(thetaAll)=="Inf")>0){
    thetaAll[iiH,"Inf"] <- 1/thetaAll[iiH,"Inf"]
  }
  if(sum(names(thetaAll)=="Vex")>0){
    thetaAll[iiH,"Vex"] <- 1/thetaAll[iiH,"Vex"]
  }
  if(sum(names(thetaAll)=="MuV")>0){
    thetaAll[iiH,"MuV"] <- 1/thetaAll[iiH,"MuV"]
  }
  
  ## Initial compartment conditions
  #theta_initAll[iiH,"fpS_init"]=239800; theta_initAll[iiH,"fpI_init"]=200;
  initial_inf=as.numeric(thetaAll[iiH,'inf0'])/2
  init_vec=as.numeric(thetaAll[iiH,'vec0']/2)
  init_rec=as.numeric(thetaAll[iiH,'rec0'])*popsizeTot
  
  theta_initAll[iiH,"r_init"]=init_rec
  theta_initAll[iiH,"e_init"]=initial_inf 
  theta_initAll[iiH,"i1_init"]=initial_inf
  theta_initAll[iiH,"em_init"]=init_vec
  theta_initAll[iiH,"im_init"]=init_vec
  
  theta_initAll[iiH,"s_init"]=popsizeTot-theta_initAll[iiH,"i1_init"]-theta_initAll[iiH,"e_init"]-theta_initAll[iiH,"r_init"]
  theta_initAll[iiH,"sm_init"]=1-theta_initAll[iiH,"em_init"]-theta_initAll[iiH,"im_init"]
  
  theta_initAll[iiH,"ed_init"]=0; theta_initAll[iiH,"id_init"]=0; theta_initAll[iiH,"t1d_init"]=0; theta_initAll[iiH,"t2d_init"]=0
  theta_initAll[iiH,"sd_init"]=(popsizeTot*(1-0.331))-theta_initAll[iiH,"id_init"]-theta_initAll[iiH,"ed_init"]-theta_initAll[iiH,"t1d_init"]-theta_initAll[iiH,"t2d_init"]
  }

## Covariance matrices 
parameters_est <- read.csv(paste0("data/",parameter_est_file,".csv"), stringsAsFactors = F) 
parms_to_est <- parameters_est$parameters_est
if(sample.start.point==T){parms_to_est <- c(parms_to_est, "t0")}
compartments_to_est <- parameters_est$compartments_est

#theta - global
nparam=length(theta) 
npc=rep(0,nparam)
npc[match(parms_to_est,names(theta))]=1
cov_matrix_theta0 = diag(npc)
  colnames(cov_matrix_theta0)=names(theta)
  rownames(cov_matrix_theta0)=names(theta)

#thetaAll - local
nparamA=length(thetaAll[1,]) 
npcA=rep(0,nparamA)
npcA[match(parms_to_est,names(thetaAll[1,]))]=1
cov_matrix_thetaAll = diag(npcA)
  colnames(cov_matrix_thetaAll)=names(thetaAll[1,])
  rownames(cov_matrix_thetaAll)=names(thetaAll[1,])

#theta_initAll 
lmv=length(theta_initAll[1,]) 
npcInit=rep(0,lmv)
npcInit[match(compartments_to_est,names(theta_initAll[1,]))]=1
cov_matrix_theta_initAll = diag(npcInit)
  colnames(cov_matrix_theta_initAll)=names(theta_initAll[1,])
  rownames(cov_matrix_theta_initAll)=names(theta_initAll[1,])

## Set up empty frames to store results
length.of.results.frames <- round(MCMC.runs/thinning.parameter,0)+1

thetatab=matrix(NA,nrow=(length.of.results.frames+1),ncol=length(theta))
colnames(thetatab)=names(theta)
thetatab[1,]=theta

thetaAlltab=array(NA, dim=c(length.of.results.frames+1,locnn,length(thetaAll[1,])),dimnames=list(NULL,NULL,names(thetaAll)))
thetaAlltab[1,,]=as.matrix(thetaAll)

theta_initAlltab=array(NA, dim=c(length.of.results.frames+1,locnn,length(theta_initAll[1,])),dimnames=list(NULL,NULL,names(theta_initAll)))
theta_initAlltab[1,,]=as.matrix(theta_initAll)

prior=rep(1,(length.of.results.frames+1))
sim_liktab=rep(-Inf,(length.of.results.frames+1))
accepttab=rep(NA,(length.of.results.frames))

full.series = seq(start.output.date,end.output.date,7)
max.length = length(full.series)

#total pop
#c
c_trace_tab=array(NA, dim=c(length.of.results.frames+1,locnn,max.length)) 
#DENV cases
cd_trace_tab=array(NA, dim=c(length.of.results.frames+1,locnn,max.length)) 
#s
s_trace_tab=array(NA, dim=c(length.of.results.frames+1,locnn,max.length)) 
#r
r_trace_tab=array(NA, dim=c(length.of.results.frames+1,locnn,max.length)) 
#x
x_trace_tab=array(NA, dim=c(length.of.results.frames+1,locnn,max.length)) 

return(list(prior=prior,
            sim_liktab=sim_liktab,
            accepttab=accepttab,
            max.length=max.length,
            c_trace_tab=c_trace_tab,
            cd_trace_tab=cd_trace_tab,
            s_trace_tab=s_trace_tab,
            r_trace_tab=r_trace_tab,
            x_trace_tab=x_trace_tab,
            thetatab=thetatab,
            thetaAlltab=thetaAlltab,
            theta_initAlltab=theta_initAlltab,
            cov_matrix_theta0=cov_matrix_theta0,
            cov_matrix_thetaAll=cov_matrix_thetaAll,
            cov_matrix_theta_initAll=cov_matrix_theta_initAll))
}
