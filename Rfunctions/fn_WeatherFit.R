#' WeatherFit likelihood fn 
#' 
#' Likelihood for seasonality parameters fitting to temperature data
#' @param amp amplitude of seasonal wave. Defaults to 0
#' @param mid midpoint of seasonal wave. Defaults to -0.048 to fix at Feb
#' @param time.vals Time series over which to fit sin curve to data
#' @param beta Average value from data to act as midpoint of sincurve
#' @param data Data values
#' @keywords Seasonality
#' @export

## likelihood function when fitting model to climate data - called in weather.fit
likelihood_fn <- function(amp, mid, time.vals=time.vals, beta=beta, data=data){
  j=1; sincurve=NULL
  for(i in time.vals){
    sincurve[j] <- beta*(1 + amp * sin(((i/365.25) - mid)*2*pi))
    j=j+1
  }
  lik <- sum(log(dnorm(sincurve, mean=data, sd=sd(data, na.rm=T))), na.rm = T)
  if(is.na(lik)){lik=0}
  if(lik==-Inf){lik=0}
  return(lik)
}

#' WeatherFit 
#' 
#' Fit seasonality parameters to temperature data
#' @param time.length timeseries to fit over
#' @param timeseries timeseries to fit over
#' @param data.series data series to fit to 
#' @param iter number of iterations to fit over
#' @param tuning Tuning param for proposal distributions
#' @keywords Seasonality
#' @export

weather.fit <- function(data, iter, tuning){
  time.vals <- seq(1,length(data),1)
  beta <- mean(data, na.rm = T)
  delta <- tuning
  
  #init
  amp.cur <- 0.5
  mid.cur <- 0.5
  
  # init priors
  amp.prior.cur <- priorBeta_amp(amp.cur)
  mid.prior.cur <- priorBeta_mid(mid.cur)
  
  # init lik
  lik.cur.amp <- likelihood_fn(amp.cur,mid.cur, time.vals=time.vals, beta=beta, data=data)
  lik.cur.mid <- likelihood_fn(amp.cur,mid.cur, time.vals=time.vals, beta=beta, data=data)

  amp <- rep(NA, iter)
  amp[1] <- amp.cur
  mid <- rep(NA, iter)
  mid[1] <- mid.cur
  accept.tab <- matrix(NA, nrow=2, ncol=iter)
  lik <- matrix(NA, nrow=2, ncol=iter)
  epsilon0.amp=0.1
  epsilon0.mid=0.1
  accept.rate.amp=0.234
  accept.rate.mid=0.234
  
  # mcmc loop for AMP
  for(i in 2:iter){
    
    epsilon0.amp = max(min(0.1,exp(log(epsilon0.amp)+(accept.rate.amp-0.234)*0.999^i)),1e-6) # Stop epsilon getting too big or small
    epsilon0.mid = max(min(0.1,exp(log(epsilon0.mid)+(accept.rate.mid-0.234)*0.999^i)),1e-6) # Stop epsilon getting too big or small
    
    #proposal P
    #amp.prop <- runif(1,0,1)
    amp.prop <- rnorm(1, mean=amp.cur, sd=epsilon0.amp)
    if(amp.prop>0){
    #likelihoods
    lik.prop.amp <- likelihood_fn(amp=amp.prop, mid=mid.cur, time.vals=time.vals, beta=beta, data=data)
    
    #priors
    amp.prior.prop <- priorBeta_amp(amp.prop)
    
    #MH number
    MH.alg <- (amp.prior.prop/amp.prior.cur)*           ## prior
              exp(lik.prop.amp-lik.cur.amp)*            ## lik
              (dunif(amp.prop,0,1)/dunif(amp.cur,0,1))  ##q ratio
    
    #selection decision
    if(runif(1) < min(1, MH.alg)){
      lik.cur.amp <- lik.prop.amp
      lik[1,i] <- lik.prop.amp
      amp.prior.cur <- amp.prior.prop
      amp[i] <- amp.prop 
      amp.cur <- amp.prop
      accept.tab[1,i] <- 1
    }else{
      amp[i] <- amp.cur
      accept.tab[1,i] <- 0
      lik[1,i] <- lik.cur.amp
    }
    }else{
      amp[i] <- amp.cur
      accept.tab[1,i] <- 0
      lik[1,i] <- lik.cur.amp
    }
    
    # repeat for MID
    mid.prop <- rnorm(1, mean=mid.cur, sd=epsilon0.mid)
    mid.prop=min(mid.prop, 2-mid.prop)
    if(mid.prop>0){
      #likelihoods
      lik.prop.mid <- likelihood_fn(amp=amp.cur, mid=mid.prop, time.vals=time.vals, beta=beta, data=data)
      
      #priors
      mid.prior.prop <- priorBeta_mid(mid.prop)
      
      #MH number
      MH.alg <- (mid.prior.prop/mid.prior.cur)*              ## prior
                exp(lik.prop.mid-lik.cur.mid)*               ## lik
                1 #(dunif(mid.prop,0,1)/dunif(mid.cur,0,1))  ##q ratio      
      #selection decision
      if(runif(1) < min(1, MH.alg)){
        lik.cur.mid <- lik.prop.mid
        lik[2,i] <- lik.prop.mid
        mid.prior.cur <- mid.prior.prop
        mid[i] <- mid.prop 
        mid.cur <- mid.prop
        accept.tab[2,i] <- 1
      }else{
        mid[i] <- mid.cur
        accept.tab[2,i] <- 0
        lik[2,i] <- lik.cur.mid
      } 
    }else{
      mid[i] <- mid.cur
      accept.tab[2,i] <- 0
      lik[2,i] <- lik.cur.mid
    }
  if(i<20){
    accept.rate.amp=0.234
    accept.rate.mid=accept.rate.amp
  }else{
    accept.rate.amp=sum(accept.tab[1,1:i],na.rm=T)/i
    accept.rate.mid=sum(accept.tab[2,1:i],na.rm=T)/i
  }
  }
  return(list(amp=amp , mid=mid, accept=accept.tab, lik=lik))
} # end Function