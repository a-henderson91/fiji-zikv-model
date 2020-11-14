# weekly: 1.98 mins
# Monthly with printing : 2.25 mins
# Monthly no printing : 2.24 mins
# 
# LSODA
# weekly: 1.3 mins
# monthly: 1.22 mins
# 
# With some fiddling:
# weekly:  1.5 mins
# monthly:  1.06 mins

dt <- 7
data <- load.data.multistart(Virus = "ZIKV", startdate = start.output.date, serology.file.name = serology.excel, init.values.file.name = init.conditions.excel, add.nulls = 0) #virusTab[iiH], dataTab[iiH])
  list2env(data,globalenv())
dt
zika_single_sim(0.22)


dt <- (7*52)/12
data <- load.data.multistart.month(Virus = "ZIKV", startdate = start.output.date, serology.file.name = serology.excel, init.values.file.name = init.conditions.excel, add.nulls = 0) #virusTab[iiH], dataTab[iiH])
  list2env(data,globalenv())
dt
time.vals
  zika_single_sim(0.22)


library(microbenchmark)
  
theta_test <- theta
init_test <- init1

weekly_fit_fn <- function(){
  dt <- 7
  data <- load.data.multistart(Virus = "ZIKV", startdate = start.output.date, serology.file.name = serology.excel, init.values.file.name = init.conditions.excel, add.nulls = 0) #virusTab[iiH], dataTab[iiH])
    list2env(data,globalenv())
  dt; time.vals
  zikv_model_ode(theta_test, init_test, time.vals.sim=time.vals)
}
monthly_fit_fn <- function(){
  dt <- (7*52)/12
  data <- load.data.multistart.month(Virus = "ZIKV", startdate = start.output.date, serology.file.name = serology.excel, init.values.file.name = init.conditions.excel, add.nulls = 0) #virusTab[iiH], dataTab[iiH])
    list2env(data,globalenv())
  dt; time.vals
  zikv_model_ode(theta_test, init_test, time.vals.sim=time.vals)
}

microbenchmark(
  zikv_model_ode(theta_test, init_test, time.vals.sim=time.vals),
  zikv_model_ode2(theta_test, init_test, time.vals.sim=time.vals)
  )
microbenchmark(
  weekly_fit_fn(),
  monthly_fit_fn()
  )
microbenchmark(
  load.data.multistart(Virus = "ZIKV", startdate = start.output.date, serology.file.name = serology.excel, init.values.file.name = init.conditions.excel, add.nulls = 0) ,
  load.data.multistart.month(Virus = "ZIKV", startdate = start.output.date, serology.file.name = serology.excel, init.values.file.name = init.conditions.excel, add.nulls = 0) 
)
