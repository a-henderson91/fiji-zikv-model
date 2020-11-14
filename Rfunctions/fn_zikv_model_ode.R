#' ODEs for SEIR-SEI with cross immunity from denv infection
#' 
#' This function solves an SEIR-SEI vector borne disease model.
#' @param theta Vector of parameters for model
#' @param init.state List of initial values
#' @param time.vals.sim List of time values to simulate over
#' @export

zikv_model_ode <- function(theta, init.state, time.vals.sim) {
  SIR_ode <- function(time, state, theta) {
    ## extract parameters from theta
    Nsize <-  theta[["npop"]]
    rho <- theta[['rho']]
    omega_d <- theta[['omega_d']]
    chi <- theta[['chi']]
    eta <- theta[['eta']]
    mu <- theta[['mu']]
      
      beta_d <- theta[['beta_d']]
      alpha_d <- theta[['alpha_d']]
      gamma_d <- theta[['gamma_d']]
      beta_h1 <- theta[['beta_h']]*
        seasonal_f(time, date0=theta[["shift_date"]],amp=theta[["beta_v_amp"]],mid=theta[["beta_v_mid"]])*
        control_f(time, base=theta[["beta_base"]], mid=theta[["beta_mid"]], width=theta[["beta_width"]]) 
      alpha_h <-  theta[["Exp"]]
      gamma <-    theta[["Inf"]]
      
      # FOI
      lambda_h <- beta_h1
      
      ## extract initial states from theta_init
      S <- state[["s_init"]]
      E <- state[["e_init"]]
      I <- state[["i_init"]]
      R <- state[["r_init"]]
      C <- state[["c_init"]] 
      Sd <- state[["sd_init"]]
      Ed <- state[["ed_init"]]
      Id <- state[["id_init"]]
      T1d <- state[["t1d_init"]]
      T2d <- state[["t2d_init"]]
      Cd <- state[["cd_init"]]
      
      ## extinction if not at least 1 infected
      Ipos = extinct(I,1) # Need at least one infective
      Idpos = extinct(Id,1) # Need at least one infective
      
      # Introduction of ZIKV infections
      intro_zikv <- intro_f(time, mid = theta[["zika_start_point"]], width = theta[["intro_width"]], base = theta[["intro_base"]]) 
      initDenv <- intro_f(time, mid = theta[["denv_start_point"]], width = 0.25, base = 160) ## fixed so that approx ~160 introduction happen on 2013-10-27
      
      # Human population
      Npop = S + E + I + R
      dS  =  eta*Npop - S*(lambda_h*I/Nsize)*Ipos - chi*Sd*(beta_d*Id/Nsize) + chi*(2*omega_d*T2d) - mu*S
      dE  =  S*(lambda_h*I/Nsize)*Ipos - alpha_h*E - mu*E
      dI  = alpha_h*E  - gamma*I + intro_zikv - mu*I
      dR  = gamma*I - rho*R - mu*R
      dC  = alpha_h*E
      
      # Denv infection and temporary immunity
      dSd = -Sd*(beta_d*Id/Nsize)*Idpos
      dEd = Sd*(beta_d*Id/Nsize)*Idpos - alpha_d*Ed 
      dId = alpha_d*Ed - gamma_d*Id + initDenv 
      dT1d = gamma_d*Id - 2*omega_d*T1d
      dT2d = 2*omega_d*T1d - 2*omega_d*T2d
      dCd = alpha_d*Ed
      
    return(list(c(dS,dE,dI,dR,dC,dSd,dEd,dId,dT1d,dT2d,dCd)))
  }
  traj <- as.data.frame(ode(init.state, time.vals.sim, SIR_ode, theta, method = "lsoda"))
  return(traj)
}