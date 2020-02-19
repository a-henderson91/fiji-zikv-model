#' Load data for model fitting 
#' 
#' This function loads surveillance cases, serology results and initial parameter values for vector-borne model fitting
#' @param add.nulls Number of time series values to add to end (default is 0)
#' @param startdate Start month of time series
#' @param Virus Name of virus for serology results. e.g. 'Zika'
#' @param cases.file.name Name of csv file with surveillance data
#' @param serology.file.name Name of csv file with serological data 
#' @param init.values.file.name Name of csv file with initial parameter values
#' @keywords import
#' @export

load.data.multistart <- function(Virus, add.nulls=0, startdate, serology.file.name, init.values.file.name){
  if(Virus=="DEN3"){locnn_iiH <- 2}else{locnn_iiH <- 1}
    cases.file.name <- dataTab[locnn_iiH]
  ## Load case data. First column must be 'date'. Column with data must be under same name as in locationtab
  timeser <- read.csv(paste0("data/",cases.file.name,".csv"), stringsAsFactors = F)  # Load ZIKA  data
    colID=(names(timeser)==locationtab[locnn_iiH])
    endwk=length(timeser[,1])
    startdate=min(as.Date(timeser$date,origin="1970-01-01"), startdate)
    enddate=max(as.Date(timeser$date,origin="1970-01-01")+add.nulls)
    pickdate=(as.Date(timeser$date,origin="1970-01-01")>=startdate & as.Date(timeser$date,origin="1970-01-01")<=enddate)
    y.vals.null <- rep(0, round((min(as.Date(timeser$date,origin="1970-01-01")) - startdate)/7))
    y.vals.add <- rep(0, (add.nulls/7))
    y.vals = c(y.vals.null, as.numeric(timeser[,colID])[pickdate], y.vals.add)
    endwk=endwk+length(y.vals.null)
    firstentry=min(c(1:endwk)[])
    y.vals=y.vals[firstentry:endwk]
    if(length(y.vals.null)==0){
      date.vals.null=NULL
    }else{
      date.vals.null=seq(startdate, startdate+((length(y.vals.null)-1)*7), 7)
    }
    if(add.nulls==0){
      date.vals.add=NULL
    }else{
      date.vals.add =seq(round(max(as.Date(timeser$date, origin="1970-01-01"))), round(enddate), 7)
    }
    ## adjust date vals by difference between new.start.time and original specified
    date.vals <- as.Date(timeser$date)
    adj <- length(y.vals.null) - (min(as.Date(timeser$date,origin="1970-01-01")) - startdate)/7
    date.vals <- date.vals + (adj*7)
    date_list=c(date.vals.null, date.vals[pickdate], date.vals.add)
    time.vals = as.numeric(date_list-min(date_list)) + (firstentry-1) * 7 # Shift by start date so both time series line up -- IMPORTANT FOR SEASONALITY
    time.vals = time.vals+time.vals[2]-time.vals[1]
    date.vals = as.Date(time.vals, origin = startdate) 
  
  ## Load serology data
  serology <- read.csv(paste0("data/",serology.file.name,".csv"), stringsAsFactors = F)  # Load ZIKA  data
    nLUM <- c(serology$sero2013[serology$virus==Virus & serology$age == 'c' & serology$serology == 'pos']+
               serology$sero2013[serology$virus==Virus & serology$age == 'a' & serology$serology == 'pos'],
                  serology$sero2015[serology$virus==Virus & serology$age == 'c' & serology$serology == 'pos']+
                  serology$sero2015[serology$virus==Virus & serology$age == 'a' & serology$serology == 'pos'],
                  serology$sero2017[serology$virus==Virus & serology$age == 'c' & serology$serology == 'pos']+
                  serology$sero2017[serology$virus==Virus & serology$age == 'a' & serology$serology == 'pos']  
                )
    nPOP <- c(serology$sero2013[serology$virus==Virus & serology$age == 'c' & serology$serology == 'total']+
              serology$sero2013[serology$virus==Virus & serology$age == 'a' & serology$serology == 'total'],
                  serology$sero2015[serology$virus==Virus & serology$age == 'c' & serology$serology == 'total']+
                  serology$sero2015[serology$virus==Virus & serology$age == 'a' & serology$serology == 'total'],
                  serology$sero2017[serology$virus==Virus & serology$age == 'c' & serology$serology == 'total']+ 
                  serology$sero2017[serology$virus==Virus & serology$age == 'a' & serology$serology == 'total']  
                )
  
  ## Load initial parameter values
  thetaR_IC_global <- read.csv(paste0("data/",init.values.file.name,"_global.csv"),stringsAsFactors=FALSE)
  thetaR_IC_local <- read.csv(paste0("data/",init.values.file.name,"_local.csv"),stringsAsFactors=FALSE)
  
  ## Set prior distribution parameters 
  return(list(y.vals=y.vals, time.vals=time.vals, date.vals=date.vals, nLUM=nLUM, nPOP=nPOP, 
              thetaR_IC_global=thetaR_IC_global, thetaR_IC_local=thetaR_IC_local))
}
