##' Compute prior
##'
##' @param iiH Numeric indicator of location given by location vector "locationtab[]"
##' @param thetaAlltab Current values of parameter set theta
##' @param thetaAllstar Proposed values of parameter set theta
##' @export
# thetaAlltab=thetaAlltab_current[iiH,];  thetaAllstar=thetaA_star; covartheta = cov_matrix_thetaA
ComputePrior <- function(iiH, thetaAlltab, thetaAllstar, covartheta){
  # Compute prior of proposed parameter set theta
  p_theta_star =  priorChi(thetaAllstar["chi"])*
                  priorEpsilon(thetaAllstar["epsilon"])*
                  priorAlpha(thetaAllstar["alpha"])*
                  priorBeta(thetaAllstar["beta_h"])*
                  priorInitInf(thetaAllstar["intro_base"])*
                  priorIntro(thetaAllstar["intro_mid"])

  # If vector control is included in the model - include prior density of proposed
  #   strength of control measure
  if(vector.control==T){
    #p_theta_star <- p_theta_star*prior_control(c(thetaAllstar["beta_base"], thetaAllstar["beta_mid"]))
    p_theta_star <- p_theta_star*prior_vectorcontrol(thetaAllstar["beta_base"])
  }
  # If seasonal transmission is estimated, include in prior
  if(seasonal.transmission==T){
    p_theta_star <- p_theta_star*
      priorBeta_amp(thetaAllstar["beta_v_amp"])*
      priorBeta_mid(thetaAllstar["beta_v_mid"])
  }
  # Compute prior of current parameter set theta
  p_theta = priorChi(thetaAlltab["chi"])*
            priorEpsilon(thetaAlltab["epsilon"])*
            priorAlpha(thetaAlltab["alpha"])*
            priorBeta(thetaAlltab["beta_h"])*
            priorInitInf(thetaAlltab["intro_base"])*
            priorIntro(thetaAlltab["intro_mid"])

  # If vector control is included in the model - include prior density of proposed
  #   strength of control measure
  if(vector.control==T){
    #p_theta <- p_theta*prior_control(c(thetaAlltab["beta_base"], thetaAlltab["beta_mid"]))
    p_theta <- p_theta*prior_vectorcontrol(thetaAlltab["beta_base"])
  }
  # If seasonal transmission is estimated, include in prior
  if(seasonal.transmission==T){
    p_theta <- p_theta*
      priorBeta_amp(thetaAlltab["beta_v_amp"])*
      priorBeta_mid(thetaAlltab["beta_v_mid"])
  }
  # formatting and replace NAs with 0
  names(p_theta_star)=NULL;names(p_theta)=NULL
  if(is.na(p_theta_star)==T){p_theta_star=1}
  if(is.na(p_theta)==T){p_theta=1}
  return(list(prior.star=p_theta_star,prior=p_theta))
}
