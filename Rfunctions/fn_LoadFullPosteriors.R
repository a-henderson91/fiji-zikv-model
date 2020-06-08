##' load full posteriors for convergence checks
##' 
##' Load and combine posterior chains following MCMC run
#' @param agestructure Binary indicator variable if model is age structured (1) or not (0) between children and adults. Defaults to NULL in which case no sampling happens 
#' @param iiH Locationtab index
#' @param mcmc.burn Between 0 and 1. Proportion of iterations to discard.
##' @export 
##' 
# load.run.name <- run.name; file.path <- "posterior_output"; iiH; mcmc.burn <- mcmc.burn

load.full.posteriors <- function(load.run.name, file.path="posterior_outputZ", iiH, mcmc.burn){
  m.tot <- length(list.files(path = paste0(file.path,"/"), pattern=paste0("*",load.run.name)))
  # iiM=1
  load_theta <- function(iiM, load.run.name){
    load(paste0(file.path,"/outputR_",iiM,"_",load.run.name,".RData",sep=""))
    theta_select <- data.frame(thetaAlltab[,iiH,])
    
    # Return a character vector of variable names which have 0 variance
    theta_select_names <- names(theta_select)[vapply(theta_select, function(x) var(x, na.rm = T) > 1e-10, logical(1))]
    
    theta_select <- theta_select[, theta_select_names]
    mcmc_samples <- length(sim_liktab)
    maxB <- sum(sim_liktab!=-Inf)/mcmc_samples
    picks <- c(1:round(maxB*mcmc_samples))
    
    list(trace = theta_select[picks,], acceptance.rate = sum(accepttab, na.rm = T)/length(accepttab))
  }
  theta_1 <- load_theta(1, load.run.name)
  theta_2 <- load_theta(2, load.run.name)
  
  picks <- c(1:length(load_theta1[,1]))
  
  return(list(picks=picks,
              thetatab=thetatab,
              sim_liktab=sim_liktab,
              accepttab=accepttab
  ))
}
