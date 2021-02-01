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

load.data.multistart.month <- function(Virus, add.nulls=0, startdate, serology.file.name, init.values.file.name){
  if(Virus=="DEN3"){locnn_iiH <- 2}else if(Virus=="ZIKV"){locnn_iiH <- 1}
    cases.file.name <- dataTab[locnn_iiH]
  ## Load case data. First column must be 'date'. Column with data must be under same name as in locationtab
  timeser <- read.csv(paste0("data/",cases.file.name,".csv"), stringsAsFactors = F)  # Load ZIKA  data
    colID=(names(timeser)==locationtab[locnn_iiH])
    endwk=length(timeser[,1])
    startdate <- min(as.Date(timeser$date,origin="1970-01-01"), startdate)
    enddate <- max(as.Date(timeser$date,origin="1970-01-01"))
    pickdate <- (as.Date(timeser$date,origin="1970-01-01")>=startdate & as.Date(timeser$date,origin="1970-01-01")<=enddate)
    y.vals = as.numeric(timeser[,colID])[pickdate]
    firstentry=min(c(1:endwk)[])
    y.vals=y.vals[firstentry:endwk]
    
    ## adjust date vals by difference between new.start.time and original specified
    date.vals <- as.Date(timeser$date)
    date_vals_prefix <- seq.Date(from = startdate, to = min(date.vals), by = "7 days")
    date_list <- c(date_vals_prefix, date.vals)
    time.vals = as.numeric(date_list-min(date_list)) + (firstentry-1) * 7 # Shift by start date so both time series line up -- IMPORTANT FOR SEASONALITY
    time.vals = time.vals+time.vals[2]-time.vals[1]
    date.vals = as.Date(time.vals, origin = startdate) 
    
    ### summarise by Month
    y.vals.full <- c(rep(0, length(date_vals_prefix)), y.vals)
    new_data <- bind_cols(
      "date.vals" = date.vals,
      "y.vals.full" = y.vals.full
      ) %>%
      mutate(month = as.numeric(lubridate::month(date.vals)),
             year = as.numeric(lubridate::year(date.vals))) %>%
      group_by(year, month) %>%
      summarise(date.vals = min(date.vals),
                y.vals = sum(y.vals.full, na.rm = T),
                .groups = "drop") 
    new_data$time.vals <- as.numeric(new_data$date.vals-min(new_data$date.vals)) + (firstentry-1) * 30 # Shift by start date so both time series line up -- IMPORTANT FOR SEASONALITY

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
  thetaR_IC_global <- read.csv(paste0("data/thetaR_IC_global.csv"),stringsAsFactors=FALSE)
  thetaR_IC_local <- read.csv(paste0("data/",init.values.file.name,"_local.csv"),stringsAsFactors=FALSE)
  
  ## Set prior distribution parameters 
  return(list(y.vals=new_data$y.vals, time.vals=new_data$time.vals, date.vals=new_data$date.vals, nLUM=nLUM, nPOP=nPOP, 
              thetaR_IC_global=thetaR_IC_global, thetaR_IC_local=thetaR_IC_local))
}
