#' calculate_r0
#'
#' @param th_in Binary indicator variable if model is age structured (1) or not (0) between children and adults. Defaults to NULL in which case no sampling happens 
#' @param sus_c Suceptible child human hosts
#' @param sus_a Suceptible adult human hosts
#' @param sm_c Suceptible child moisquitoes
#' @param sm_a Suceptible adult moisquitoes
#' @param b_vary Defaults to 1
#' @export

calculate_r0 <- function(th_in,sus_c=1,sus_a=1,sm_c=1,sm_a=1,b_vary=1,control=1){
  # Rate humans get infected -- FORMULATION WITH ONE MOSQUITO POP AND DIFFERENT BITING RATES
  b_hv =  b_vary * th_in$beta_h
  b_hh = 0
  
  # Rate vectors get infected
  b_vh = 0
  b_vv = 0
  
  rr_hh <- rep(0,length(b_vary)); 
  rr_vv <- rr_hh;
  
  exp_h <- th_in$Exp
  inf_p <- th_in$Inf.
  
  rr_hh = (sus_c*b_hv/(exp_h + inf_p))
  rr_vv = rr_hh
  rr_hv = rr_hh
  rr_vc = rr_hh
  
  r0_hh = b_hv/(exp_h + inf_p)
  r0_vv = rr_hh
  r0_hv = rr_hh
  r0_vc = rr_hh
  
  rr_post = NULL
  r0_post = NULL
  for(ii in 1:length(b_vary)){
    rr_post=c(rr_post,( max(Re(eigen(matrix(c(rr_hh[ii],rr_hv[ii],rr_vc[ii],rr_vv[ii]),nrow=2))$values)))  )
    r0_post=c(r0_post,( max(Re(eigen(matrix(c(r0_hh[ii],r0_hv[ii],r0_vc[ii],r0_vv[ii]),nrow=2))$values)))  )
  }
  return( list(r0_out=r0_post, rr_out=rr_post,rr_mat = matrix(c(rr_hh[1],rr_hv[1],rr_vc[1],rr_vv[1]),nrow=2,byrow=T)) )
}
