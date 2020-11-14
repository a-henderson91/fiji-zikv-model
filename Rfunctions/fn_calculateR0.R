#' calculate_r0
#'
#' @param th_in parameter set
#' @param sus Suceptible adults 
#' @param b_vary Defaults to 1
#' @export

calculate_r0 <- function(th_in,sus = 1, b_vary = 1){
  # Rate humans get infected -- FORMULATION WITH for SEIR model with demographics
  b_h =  b_vary * th_in$beta_h
   
  rr_hh <- rep(0,length(b_vary)); 
  
  exp_h <- th_in$Exp
  inf_p <- th_in$Inf.
  death_rate <- 1/(th_in$mu*365.25)
  
  r0_post = (b_h*exp_h)/((death_rate + inf_p)*(death_rate + exp_h))
  rr_post = (sus*b_h*exp_h)/((death_rate + inf_p)*(death_rate + exp_h))
  
  return( list(r0_out=r0_post, rr_out=rr_post))
}
