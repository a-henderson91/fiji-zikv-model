#' Sample the vector of parameters from a multivariate normal distribution
#'
#' Sample a new vector of parameters (either global or local) and initial conditions using a multivariate normal dist

#theta_in=thetaAlltab_current[iiH,]; theta_init_in=theta_initAlltab_current[iiH,]; covartheta=0*cov_matrix_thetaA; covartheta_init=0*cov_matrix_theta_init; global=0
#theta_in=thetaAlltab_current[iiH,]; theta_init_in=theta_initAlltab_current[iiH,]; covartheta=cov_matrix_thetaA; covartheta_init=cov_matrix_theta_init; global=0
SampleTheta<-function(theta_in, theta_init_in, covartheta, covartheta_init, global=NULL){
  ## Parameters
    # sample new parameters from nearby using multivariate normal distribution:
      mean_vector_theta = theta_in
      mean_vector_theta0 = mean_vector_theta
      theta_star = as.numeric(exp(rmvnorm(1,log(mean_vector_theta0), covartheta)))
      names(theta_star)=names(theta_in)

      if(sum(names(theta_star)=="rep")>0){ # check theta contains this vector
        theta_star[["rep"]]=min(theta_star[["rep"]],2-theta_star[["rep"]]) # Ensure reporting between zero and 1
      }
      if(sum(names(theta_star)=="chi")>0){ # check theta contains this vector
        theta_star[["chi"]]=min(theta_star[["chi"]],2-theta_star[["chi"]]) # Ensure reporting between zero and 1
      }
      if(sum(names(theta_star)=="epsilon")>0){ # check theta contains this vector
        theta_star[["epsilon"]]=min(theta_star[["epsilon"]],2-theta_star[["epsilon"]]) # Ensure reporting between zero and 1
      }
      if(sum(names(theta_star)=="inf0")>0){ # check theta contains this vector
        theta_star[["inf0"]]=max(0,min(theta_star[["inf0"]],1000)) # Ensure initial infectious < total pop
      }
      if(sum(names(theta_star)=="beta_base")>0){
        theta_star[["beta_base"]]=min(theta_star[["beta_base"]],2-theta_star[["beta_base"]]) # Ensure amplitude between zero and 1
      }
      if(sum(names(theta_star)=="intro_mid")>0){
        theta_star[["intro_mid"]]=max(0, min(theta_star[["intro_mid"]],800))
      }

  ## Initial conditions
  theta_init_star = theta_init_in

  popsizeTot=theta_init_star["s_init"]+theta_init_star["e_init"]+theta_init_star["i1_init"]+theta_init_star["r_init"]

  initial_inf=as.numeric(theta_star['inf0'])/2
  init_rec=popsizeTot*as.numeric(theta_star['rec0'])

  theta_init_star["r_init"]=init_rec
  theta_init_star["e_init"]=initial_inf; theta_init_star["i1_init"]=initial_inf

  theta_init_star["s_init"]=popsizeTot-theta_init_star["i1_init"]-theta_init_star["e_init"]-theta_init_star["r_init"]

  theta_init_star["ed_init"]=0; theta_init_star["id_init"]=0; theta_init_star["t1d_init"]=0; theta_init_star["t2d_init"]=0
  theta_init_star["sd_init"]=(popsizeTot*(1-0.331))-theta_init_star["id_init"]-theta_init_star["ed_init"]-theta_init_star["t1d_init"]-theta_init_star["t2d_init"]

  return(list(thetaS=theta_star,theta_initS=theta_init_star))
}
