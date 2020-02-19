#' Seasonal transmission rate function
#' 
#' Calculates the seasonal tranmsission rate
#' @param time t in time.vals of calc
#' @param date0 reference index value for start of transmission. Defaults to 0
#' @param amp amplitude of seasonal wave. Defaults to 0
#' @param mid midpoint of seasonal wave. Defaults to -0.048 to fix at Feb
#' @keywords Seasonality
#' @export
#' @examples 
#' seasonal_f(date0=0, amp=theta[["beta_v_amp"]], mid=theta[["beta_v_mid"]])

seasonal_f <- function(time, date0=0, amp=0, mid=pi*(3/4)){
  (1 + amp *sin(((time/365.25) - mid)*2*pi))
}

#' Flexible decline function transmission rate function
#' 
#' Sigmoid function to allow drop in transmission rate
#' @param time t in time.vals of calc
#' @param date0 reference index value for start of transmission. Defaults to 0
#' @param mask Strength of reduction in transmission rate. Defaults to 0 (i.e. no effect)
#' @param base Height of starting point of function. Defaults to 1
#' @param grad Strength of decline in transmission once begun. Defaults to 0.5
#' @param mid Point at which decline can begin. Defaults to 0
#' @keywords Sigmoid Decline
#' @export
#' @examples 

intro_f <- function(time,mid,width,base){
  xx <- 1 - (4*base)*exp(-(time-mid)/width)/(1+exp(-(time-mid)/width))^2
  if(is.nan(xx)){xx <- 1}
  xx
}

#' Flexible decline function transmission rate function
#' 
#' Sigmoid function to allow drop in transmission rate
#' @param time t in time.vals of calc
#' @param date0 reference index value for start of transmission. Defaults to 0
#' @param mask Strength of reduction in transmission rate. Defaults to 0 (i.e. no effect)
#' @param base Height of starting point of function. Defaults to 1
#' @param grad Strength of decline in transmission once begun. Defaults to 0.5
#' @param mid Point at which decline can begin. Defaults to 0
#' @keywords Sigmoid Decline
#' @export
#' @examples 
#' 
control_f <- function(time,mask=0,base=1,grad=.1,mid=0, mid2=1, width=30){
  c1 <- 1 - (4*base)*exp(-(time-mid)/width)/(1+exp(-(time-mid)/width))^2
  c1 <- c1*(time<mid)
  c2 <- 1-(base)*(1-1/(1+exp(-grad*(time-mid2))))
  c2 <- c2*(time>=mid)
  c1+c2
}

#' Reporting cases function
#' 
#' Given a number of infected hosts, this function estimates the resulting number of cases
#' @param cases List of cases from model
#' @param rep Reporting rate used to generate those cases
#' @param repvol Overdispersion parameter used in Likelihood to estimate those cases
#' @export

ReportC<-function(cases, rep, repvol){
  mu00=cases
  mu01=sapply(mu00,function(x){max(x,0)})
  sapply(mu01,function(x){rnbinom(1, mu=rep*x,size=1/repvol)})
}

#' Reporting cases function
#' 
#' Given a number of infected hosts, this function estimates the resulting number of cases
#' @param x number to test
#' @param x0 value to be greater than
#' @export

extinct <- function(x,x0){as.numeric(x>=x0)} 
