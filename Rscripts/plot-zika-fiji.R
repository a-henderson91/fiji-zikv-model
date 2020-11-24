# - - - - - - - - - - - - - - - - - - - - - - - 
# Plots and analysis of posteriors for Fiji Zika transmission model  
# Author: Alasdair Henderson
# github.com/a-henderson91/fiji-zikv-model
# - - - - - - - - - - - - - - - - - - - - - - - 

virus <- "DEN3"  # DEN3 or ZIKV

output_simulations <- F
output_diagnostics <- F
ll <- 1

here::here()
if(virus=="ZIKV"){
  iiH <- 1
  preamblecode <- "preamble_zikvfiji.R"
  post_file_path <- "posterior_output"
}else if(virus=="DEN3"){
  iiH <- 2
  preamblecode <- "preamble_denv3.R"
  post_file_path <- "posterior_denv2014fit"
}
source(here::here("Rscripts", preamblecode))
if(virus == "DEN3"){model1_name <- run.name}

# load vectorbornefit code ------------------------------------------------
all_files <- list.files(here::here("Rfunctions"))
sapply(all_files, function(xx){source(here::here("Rfunctions/", xx))})

# import posteriors -------------------------------------------------------
labelN <- 1
m.tot <- length(list.files(path = here::here(post_file_path), pattern = paste0("*", run.name)))
load_posterior_1 <- load.posteriors(load.run.name = run.name, file.path = post_file_path, iiH, mcmc.burn = mcmc.burn)
  list2env(load_posterior_1,globalenv())
  
load_DENVfit_1 <- load.posteriors(load.run.name = model1_name, file.path = "posterior_denv2014fit", iiH = 2, mcmc.burn = mcmc.burn)

# load data --------------------------------------------------------------
data <- load.data.multistart(Virus = virusTab[iiH], startdate = start.output.date, serology.file.name = serology.excel, init.values.file.name = init.conditions.excel, add.nulls = 0) #virusTab[iiH], dataTab[iiH])
  list2env(data,globalenv())
  
# load BEAST tmrca posterior - from "Export start time Central Division.R"
load("beast/fj-c-tmrca.RData")
tmrcaBEAST <- fj.C.tmrca.date

# Fig 1A - introduction time histogram ------------------------------------
tmrca <- cbind.data.frame(startdate + thetatab$intro_mid, c(startdate + tmrcaBEAST[1:length(thetatab$intro_mid)]))
  names(tmrca) <- c("posterior","prior")
  tmrca <- gather(tmrca, data_type, density)

ggplot(data=tmrca, aes(colour=data_type, fill=data_type)) +
  geom_histogram(aes(density), bins = 60) +
  labs(x = "Virus introduction date", y = "Density", title = "") +
  #scale_x_date(limits = c(startdate, as.Date("2015-12-31")), date_breaks = "3 month", date_labels = "%m/%y") +
  scale_colour_manual(labels = c("Posterior", "Prior"), values=c(col1, col2)) +
  scale_fill_manual(labels = c("Posterior", "Prior"), values=c(col1a, col2a)) +
  theme_zika_fiji() +
  theme(legend.title=element_blank())
ll <- ll+1
dev.copy(pdf, paste0("output/fig1A_introHist_",virus,".pdf"), 8, 6)
  dev.off()

# fig1B - case data -------------------------------------------------------
#virusTab[1] <- virusTab[2]
#locationtab[1] <- locationtab[2]
#dataTab[1] <- dataTab[2]
### Load DENV data and store in objects
#MCMC.runs <- min(MCMC.runs, 10000)
#  #denv <- load.data.multistart(startdate=as.Date("2013-10-27"), add.nulls=0,  "DEN3", "Central_2014_timeseries", "Fiji_serology_APR19", "thetaR_IC_zika_denvlike")
#
denv_data <- load.data.multistart(Virus = "DEN3", startdate = start.output.date, serology.file.name = serology.excel, init.values.file.name = init.conditions.excel, add.nulls = 0) #virusTab[iiH], dataTab[iiH])
  denv.cases <- denv_data$y.vals
  denv.dates <- denv_data$date.vals
## load denv posteriors
#load("data_sets/theta_denv_posteriors.RData")
#d_trace_tab <- DENV.mcmc.run$c_trace_tab[,1,]
#d_thetatab <- DENV.mcmc.run$thetaAlltab[,1,]

## build dengue and zika case dataset
data.denv <- cbind.data.frame(dates=denv.dates, denv.cases)
data.zika <- cbind.data.frame(dates=date.vals, y.vals)
DFdenv <- data.table(data.denv)
DFzikv <- data.table(data.zika)
case.data <- DFdenv[DFzikv, on=.(dates = dates), roll=-Inf] ## match on nearest date
case.data$y.vals[case.data$y.vals==0] <- NA
case.data$denv.cases[case.data$denv.cases==0] <- NA
  
## cases collected
casecollect <- read.csv(here::here("data/SamplesCollected.csv"), stringsAsFactors = F)
datesample <- as.Date(casecollect$DATINTRVW)
datesample <- datesample[datesample>as.Date("2013-10-01")]
monthsample <- (format(as.Date(casecollect$DATINTRVW), "%Y-%m"))
monthsample <- data.frame(mth=monthsample)
monthsample <- monthsample %>%
  group_by(mth) %>%
  summarize(bymth = n())
monthsample$mth <- as.character(monthsample$mth)
monthsample$mth <- as.Date(paste0(monthsample$mth,"-01"), format="%Y-%m-%d")
casecollect$x <- 1
casecollect <- casecollect %>% group_by(DATINTRVW) %>% summarise(dailysample=sum(x))
casecollect$DATINTRVW <- as.Date(casecollect$DATINTRVW) 

cexplot <- 2.5
par(mfrow=c(1,1), mar = c(3,4,1,5))
if(virus=="ZIKV"){tt="h"}else{tt="l"}
plot(casecollect$DATINTRVW, casecollect$dailysample, col="gray", lwd=cexplot, pch=16, bty='n',
     ylim=c(0,50), type="h", 
     xlab="Date", ylab="", yaxt='n', main="B", adj=0, xaxt='n')
par(new=T)
hist(c(startdate + tmrcaBEAST[1:length(thetatab$intro_mid)]), breaks = 1000, col=col2a, border=col2a,
     ylim=c(0,0.01), main="", ylab="", yaxt='n', xaxt="n", xlab="", 
     freq=F,
     xlim=c(as.Date("2013-01-08"), as.Date("2018-01-01")))
par(new=T)
if(virus == "ZIKV"){
  plot(case.data$dates, case.data$y.vals, col=col1, lwd=cexplot, pch=16, bty='n',ylim=c(0,10),type=tt, 
       xlab="Date", ylab="", yaxt='n', main="B", adj=0, xaxt='n')
  text(as.Date("2015-01-07"), 4, "tMRCA", col = col2, font = 1)
  mtext("Zika cases", side=4, line=2, col=col1, cex=1, font=1)
}
axis(4, bty='l', col.ticks = 1, col=1, col.axis=1, col.lab=1)
mtext("Serological samples collected", side=4, line=2, col="Gray", cex=1, padj = 1.5, font=2)
par(new=T)
plot(case.data$dates, case.data$denv.cases, col=col4, lwd=cexplot, pch=16,type='l',
     #xaxt="n",
     xlab="Date",yaxt="n",ylab="",frame.plot=F)
axis(2, bty='l', col.ticks = 1, col=1, col.axis=1, col.lab=1)
mtext("Dengue cases", side=2, line=2, col=col4, cex=1, font=2)
grid(NA,NULL, lty = 1, col = colgrid) 
mtext("C",side=3, adj=0, font=1)
  ll <- ll+1
dev.copy(pdf, paste0("output/fig1B_caseData_",virus,".pdf"), width = 6, height = 4)
  dev.off()

# Fig 2 model trajectories ---------------------------------------------
tMax <- dim(c_trace_tab)[2]
btsp <- 1000
svector <- matrix(NA,nrow=btsp,ncol=tMax)
cvector <- matrix(NA,nrow=btsp,ncol=tMax)
ivector <- matrix(NA,nrow=btsp,ncol=tMax)
rvector <- matrix(NA,nrow=btsp,ncol=tMax)
introvector <- matrix(NA,nrow=btsp,ncol=tMax)
cvectorDENV <- matrix(NA,nrow=btsp,ncol=tMax)
ln_denv <- length(denv.cases)
first_zikv <- min(which(y.vals>0))

set.seed(0517365978236587236)
for(ii in 1:btsp){
  pick <- sample(picks, 1)
  if(!is.na(thetatab[pick,'epsilon'])){
    epsilon <- thetatab[pick,'epsilon']
  }else{
    epsilon=0
  }
  cvector[ii,] <- ReportC(c_trace_tab[pick,1:tMax],thetatab[pick,'rep'],thetatab[pick,'repvol'])
  svector[ii,] <- s_trace_tab[pick,1:tMax]
  ivector[ii,] <- c_trace_tab[pick,1:tMax]
  rvector[ii,] <- (r_trace_tab[pick,1:tMax]/thetatab[pick,]$npop) + ((1- (r_trace_tab[pick,1:tMax]/thetatab[pick,]$npop))*epsilon)
  introvector[ii,] <- sapply(time.vals[1:tMax], function(xx){intro_f(xx, mid = thetatab[pick,'intro_mid'], width = thetatab[pick,'intro_width'], base = thetatab[pick,'intro_base'])})
}

## plot posterior of cases
tMaxDenv <- length(load_DENVfit_1$c_trace_tab[1,])
denvvector <- matrix(NA,nrow=btsp,ncol=tMaxDenv)
for(ii in 1:btsp){
  pick <- sample(load_DENVfit_1$picks, 1)
  denvvector[ii,] <- ReportC(load_DENVfit_1$c_trace_tab[pick,1:tMaxDenv],load_DENVfit_1$thetatab[pick,'rep'],load_DENVfit_1$thetatab[pick,'repvol'])
}
medD <- apply(denvvector,2,function(x){median(x, na.rm=T)})
ciD1 <- apply(denvvector,2,function(x){quantile(x,0.025, na.rm=T)})
ciD2 <- apply(denvvector,2,function(x){quantile(x,0.975, na.rm=T)})

# Estimated number of cases 
medP <- apply(cvector,2,function(x){median(x, na.rm=T)})
ciP1 <- apply(cvector,2,function(x){quantile(x,0.025, na.rm=T)})
ciP2 <- apply(cvector,2,function(x){quantile(x,0.975, na.rm=T)})
ciP150 <- apply(cvector,2,function(x){quantile(x,0.25, na.rm=T)})
ciP250 <- apply(cvector,2,function(x){quantile(x,0.75, na.rm=T)})

# Proportion recovered
medP_R <- apply(rvector,2,function(x){median(x, na.rm=T)}); 
ciP1_R <- apply(rvector,2,function(x){quantile(x,0.025, na.rm=T)}); 
ciP2_R <- apply(rvector,2,function(x){quantile(x,0.975, na.rm=T)}); 
ciR150 <- apply(rvector,2,function(x){quantile(x,0.25, na.rm=T)})
ciR250 <- apply(rvector,2,function(x){quantile(x,0.75, na.rm=T)})

# Number of infected
med_Inf <- apply(ivector,2,function(x){median(x, na.rm=T)})
ci_inf1 <- apply(ivector,2,function(x){quantile(x,0.025, na.rm=T)})
ci_inf2 <- apply(ivector,2,function(x){quantile(x,0.975, na.rm=T)})
ci_inf150 <- apply(ivector,2,function(x){quantile(x,0.25, na.rm=T)})
ci_inf250 <- apply(ivector,2,function(x){quantile(x,0.75, na.rm=T)})

# susceptible dynamics
medS <- apply(svector,2,function(x){median(x, na.rm=T)})
ciS1 <- apply(svector,2,function(x){quantile(x,0.025, na.rm=T)})
ciS2 <- apply(svector,2,function(x){quantile(x,0.975, na.rm=T)})

# intro dynamics
med_intro <- apply(introvector,2,function(x){median(x, na.rm=T)})
ci_intro1 <- apply(introvector,2,function(x){quantile(x,0.025, na.rm=T)})
ci_intro2 <- apply(introvector,2,function(x){quantile(x,0.975, na.rm=T)})

# Fig 2A - infections and cases ---------------------------------------------
par(mfrow=c(2,2))
par(mgp=c(2,0.7,0),mar = c(3,4,1.2,4))
#
y.vals.plot <- y.vals[1:tMax]
y.vals.plot[y.vals.plot==0] <- NA
dataframe.p1 <- data.table(date.vals=date.vals, y.vals.plot, 
                           medP, ciP1, ciP2, ciP150, ciP250, 
                           medP_R, ciP1_R, ciP2_R, ciR150, ciR250,
                           med_intro, ci_intro1, ci_intro2,
                           medD = c(medD, rep(NA, tMax+1-tMaxDenv)), 
                           ciD1 = c(ciD1, rep(0, tMax+1-tMaxDenv)), 
                           ciD2 = c(ciD2, rep(0, tMax+1-tMaxDenv)))
plot_fig2A <- DFdenv[dataframe.p1, on=.(dates = date.vals), roll=-Inf]

if(virus == "ZIKV"){tt = "h"}else{tt = "l"}
ylims <- c(0, max(plot_fig2A$ciP2)*1.1)
plot(plot_fig2A$dates, plot_fig2A$y.vals.plot, col=col1, cex=0.8, pch=16, bty='n',
     ylim=ylims, type = tt, lwd=2,
     xlab="Date", ylab="", yaxt='n',
     xaxt="n")
axis.Date(1, at=seq(min(plot_fig2A$dates), max(plot_fig2A$dates), by="years"), 
          labels=format(seq(min(plot_fig2A$dates), max(plot_fig2A$dates), by="years"),"%Y"), las=1)
lines(plot_fig2A$dates,plot_fig2A$medP, col=col1, lty=2)
polygon(c(plot_fig2A$dates, rev(plot_fig2A$dates)),
        c(plot_fig2A$ciP1, rev(plot_fig2A$ciP2)), col=col1a, lty=0)

  if(virus == "DEN3"){
    axis(side=2,  bty='l', col.ticks = col1a, col=col1, col.axis=col1, col.lab=col1a)
    mtext("Dengue-3 cases", side=2, line=2, col=col1, cex=1)
    grid(NA,NULL, lty = 1, col = colgrid) 
  }else if(virus == "ZIKV"){
    axis(side=4,  bty='l', col.ticks = 1, col=1, col.axis=1, col.lab=1)
    mtext("Zika cases", side=4, line=2, col=col1)
    par(new=T)
    plot(plot_fig2A$dates, plot_fig2A$ciD2, col=0, cex=0.8, pch=16, type='l',
         xaxt="n",yaxt="n",xlab="",ylab="",frame.plot=F)
    lines(plot_fig2A$dates, plot_fig2A$denv.cases, col=col2, cex=0.8, lwd = 2)
    
    lines(plot_fig2A$dates, plot_fig2A$medD, col=col2, cex=0.8, lwd = 2, lty = 2)
    polygon(c(plot_fig2A$dates, rev(plot_fig2A$dates)),
            c(plot_fig2A$ciD1, rev(plot_fig2A$ciD2)), lty=0, col=col2a)
    mtext(LETTERS[1],side=3, adj=0, font=2)
    
    axis(2, bty='l', col.ticks = 1, col=1, col.axis=1, col.lab=1)
    mtext("Dengue-3 cases", side=2, line=2, col=col2) # Label for 2nd axis
      grid(NA,NULL, lty = 1, col = colgrid) 
  }
mtext("A",side=3, adj=0, font=2)

# Fig 2B - serology and intro dynamics ---------------------------------------------
i=1; lci=NULL;uci=NULL;lci_A=NULL;uci_A=NULL;points=NULL;points_A=NULL; date=NULL
  for(date in seroposdates){
    binomtest <- binom.test(x=nLUM[i],n=nPOP[i])
    lci[i] <- binomtest$conf.int[1]
    uci[i] <- binomtest$conf.int[2]
    points[i] <- nLUM[i]/nPOP[i]
    i <- i+1
  }
dataframe.sero <- data.frame(points,lci,uci, seroposdates)
dataframe.p2 <- data.frame(date.vals=date.vals[1:tMax], y.vals.plot, 
                           medP_R, ciP1_R, ciP2_R,
                           med_intro, ci_intro1, ci_intro2)

plot(dataframe.p2$date.vals, dataframe.p2$ci_intro2, col=0, ylab="", yaxt="n", xlab="Date", bty="n", xaxt="n")
axis(side=2, bty='l', col.ticks = 1, col=1, col.axis=1, col.lab=1)
axis.Date(1, at=seq(min(plot_fig2A$dates), max(plot_fig2A$dates), by="3 months"), 
          labels=format(seq(min(plot_fig2A$dates), max(plot_fig2A$dates), by="3 months"),"%b-%y"), las=1)
par(new=T)
lines(dataframe.p2$date.vals, dataframe.p2$med_intro, col=col7, lwd = 2)
polygon(c(dataframe.p2$date.vals, rev(dataframe.p2$date.vals)), 
        c(dataframe.p2$ci_intro1,rev(dataframe.p2$ci_intro2)), lty=0, col=col7a)
mtext("Introductions", side=2, line=2, col=col7)
par(new=T)
ylims <- c(0, max(dataframe.p2$ciP2_R)*1.1)
plot(dataframe.p2$date.vals, dataframe.p2$medP_R, type='l', lty=2, col=col5a, ylim=ylims, lwd = 2,
     xaxt='n', xlab="", yaxt="n", ylab="", bty="n")
axis(side=4,  bty='l', col.ticks = 1, col=1, col.axis=1, col.lab=1)
mtext("Prop. seropostive", side=4, line=2, col=col5)
polygon(c(dataframe.p2$date.vals, rev(dataframe.p2$date.vals)), 
        c(dataframe.p2$ciP1_R, rev(dataframe.p2$ciP2_R)),
        col=col5a,lty=0)
points(dataframe.sero$seroposdates, dataframe.sero$points, col=col5, cex=2, pch=16)
lines(c(dataframe.sero$seroposdates[1],dataframe.sero$seroposdates[1]),
      c(dataframe.sero$lci[1],dataframe.sero$uci[1]), col=col5, lwd = 2)
lines(c(dataframe.sero$seroposdates[2],dataframe.sero$seroposdates[2]),
      c(dataframe.sero$lci[2],dataframe.sero$uci[2]), col=col5, lwd = 2)
lines(c(dataframe.sero$seroposdates[3],dataframe.sero$seroposdates[3]),
      c(dataframe.sero$lci[3],dataframe.sero$uci[3]), col=col5, lwd = 2)
grid(NA,NULL, lty = 1, col = colgrid) 
mtext("B",side=3, adj=0, font=2)

# Fig 2C - S and I dynamics ---------------------------------------------
data.frame2c <- tibble(
  "date.vals" = date.vals[1:tMax], 
  medS, ciS1, ciS2,
  med_Inf, ci_inf1, ci_inf2
  ) %>%
  mutate_at(c("med_Inf", "ci_inf1", "ci_inf2"), ~log(.)) %>%
  mutate_at(c("ci_inf1"), ~ifelse(ci_inf2 == -Inf, 0, .)) %>%
  mutate_at(c("ci_inf2"), ~ifelse(ci_inf1 == -Inf, 0, .)) %>%
  mutate_at(c("ci_inf1", "ci_inf2"), ~ifelse(. == -Inf, 0, .)) %>%
  mutate_at("med_Inf", ~ifelse(ci_inf1 == 0, NA, .)) %>%
  mutate_at("med_Inf", ~ifelse(ci_inf2 == 0, NA, .)) 

plot(data.frame2c$date.vals, data.frame2c$ciS2, col=0, ylab="", yaxt="n", xlab="Date", bty="n", xaxt="n",
     ylim = c(0, max(data.frame2c$ciS2)))
axis.Date(1, at=seq(min(plot_fig2A$dates), max(plot_fig2A$dates), by="3 months"), 
          labels=format(seq(min(plot_fig2A$dates), max(plot_fig2A$dates), by="3 months"),"%b-%y"), las=1)
lines(data.frame2c$date.vals, data.frame2c$medS, col=col4, lwd = 2)
polygon(c(data.frame2c$date.vals, rev(data.frame2c$date.vals)), 
        c(data.frame2c$ciS1,rev(data.frame2c$ciS2)), lty=0, col=col4a)
axis(side=2, bty='l', col.ticks = 1, col=1, col.axis=1, col.lab=1)
mtext("Number susceptible", side=2, line=2, col=col4, adj = 0.8 )
par(new = T)
plot(data.frame2c$date.vals, data.frame2c$med_Inf,
     col=0, lwd = 2, type = "l", bty = "n", xaxt = "n", yaxt = "n", ylab = "", xlab = "", ylim = c(0, max(data.frame2c$ci_inf2)))
lines(data.frame2c$date.vals, data.frame2c$med_Inf, col=col1, lwd = 2)
polygon(c(data.frame2c$date.vals, rev(data.frame2c$date.vals)), 
        c(data.frame2c$ci_inf1,rev(data.frame2c$ci_inf2)), lty=0, col=col1a)
aty <- log(c(1,10,100,1000,10000))
labels <- sapply(aty,function(i)
  as.expression(round(exp(i)))
)
axis(side=4, bty='l', col.ticks = 1, col=1, col.axis=1, col.lab=1, at=aty,labels=labels, las = 2)
mtext("Number infected", side=4, line=2, col=col1, adj = 0.2)
abline(h = 0, lty = 2, col = "gray70")
grid(NA,NULL, lty = 1, col = colgrid) 
mtext("C",side=3, adj=0, font=2)

# Fig 2D - Reproduction numbers ---------------------------------------------
plotCosRR <- NULL
plotCosR0 <- NULL
plotDecline <- NULL
btstrap <- sample(picks, 400, replace=T)
for(ii in 1:length(btstrap)){
  b_ii <- btstrap[ii]
  t.start <- 0
  time.V <- (1:tMax)*7
  date0 <- (startdate-date.vals[1]) %>% as.numeric() 
  
  beta_ii <- seasonal_f(time.V[1:tMax], date0, amp=thetatab[b_ii,'beta_v_amp'], mid=thetatab[b_ii,'beta_v_mid'])
  decline_ii <- control_f(time.V[1:tMax], base = thetatab[b_ii, "beta_base"], mid = thetatab[b_ii,"beta_mid"], width = thetatab[b_ii,"beta_width"])
  b_vary <- beta_ii
  
  s_pick <- s_trace_tab[b_ii,1:tMax]/thetatab$npop[b_ii] 
  
  output_rr <- calculate_r0(th_in=thetatab[b_ii,], sus = s_pick, b_vary=b_vary*decline_ii)
  output_rr_nocontrol <- calculate_r0(th_in=thetatab[b_ii,],sus = s_pick, b_vary = b_vary)
  
  start.rr <- output_rr$rr_out[min(which(output_rr$rr_out>0))]
  output_rr$rr_out[1:(min(which(output_rr$rr_out>0))-1)] <- start.rr
  output_rr$r0_out[1:(min(which(output_rr_nocontrol$rr_out>0))-1)] <- start.rr
  
  rr_post = output_rr$rr_out;  rr_post[length(rr_post)] = rr_post[length(rr_post)-1]
  r0_post = output_rr_nocontrol$r0_out; r0_post[length(r0_post)] = r0_post[length(r0_post)-1]
  decline_post= decline_ii  
  
  plotCosRR=rbind(plotCosRR,  rr_post)
  plotCosR0=rbind(plotCosR0,  r0_post)
  plotDecline=rbind(plotDecline, decline_post)
}

c.nume<-function(x){
  bp1=c(median(x),quantile(x,0.025),quantile(x,0.975))
  as.numeric(bp1)}
plotCosMRR <- apply(plotCosRR,2,c.nume) 
plotCosMR0 <- apply(plotCosR0,2,c.nume) 
plotDecline <- apply(plotDecline,2,c.nume) 

dataframe.p3 <- data.frame(date.vals=date.vals[1:tMax], 
                           medRR=plotCosMRR[1,], lciRR=plotCosMRR[2,], uciRR=plotCosMRR[3,],
                           medR0=plotCosMR0[1,], lciR0=plotCosMR0[2,], uciR0=plotCosMR0[3,],
                           medD=plotDecline[1,], lciD=plotDecline[2,], uciD=plotDecline[3,])
weather.data <- read.csv(here::here("data/suva-temp.csv"), stringsAsFactors = T) 
weather.data$lsd <- as.Date(weather.data$lsd)
data.series.max <- weather.data$max_air_temp[as.Date(weather.data$lsd)>=min(as.Date(date.vals)) & as.Date(weather.data$lsd)<=max(as.Date(date.vals))]
data.series.min <- weather.data$min_air_temp[as.Date(weather.data$lsd)>=min(as.Date(date.vals)) & as.Date(weather.data$lsd)<=max(as.Date(date.vals))]
data.series <- ((data.series.max-data.series.min)/2)+data.series.min
date.series <- weather.data$lsd[as.Date(weather.data$lsd)>=min(as.Date(date.vals)) & as.Date(weather.data$lsd)<=max(as.Date(date.vals))]
seas.weatherdata <- data.frame(date.vals=date.series, temp=data.series)
dataframe.p3 <- left_join(dataframe.p3, seas.weatherdata, by = "date.vals")

yht <- max(dataframe.p3$uciR0,na.rm=T)*1.2
plot(dataframe.p3$date.vals, dataframe.p3$medR0, type='l', col=col2, bty="n", ylim=c(0,yht),
     xlab="Date", ylab="", yaxt="n")
polygon(c(dataframe.p3$date.vals,rev(dataframe.p3$date.vals)),
        c(dataframe.p3$lciR0,rev(dataframe.p3$uciR0)),lty=0,col=col2a)
axis(side=2, bty='l', col.ticks = 1, col=1, col.axis=1, col.lab=1)
mtext("Reproduction number", side=2, col=1, line=2, cex=1)
lines(dataframe.p3$date.vals, dataframe.p3$medRR, col=col7)
polygon(c(dataframe.p3$date.vals,rev(dataframe.p3$date.vals)),
        c(dataframe.p3$lciRR,rev(dataframe.p3$uciRR)),lty=0,col=col7a)
abline(h=1, lty=2)    
#lines(dataframe.p3$date.vals, dataframe.p3$medD, col = col3)
#polygon(c(dataframe.p3$date.vals,rev(dataframe.p3$date.vals)),
#        c(dataframe.p3$lciD,rev(dataframe.p3$uciD)),lty=0,col=col3a)
par(new=T)
plot(dataframe.p3$date.vals, dataframe.p3$temp, type='l', col=rgb(0,0,0.1,alpha=0.4), bty="n", 
     xlab="Date", ylab="", yaxt="n")
axis(side=4, bty='l', col.ticks = 1, col=1, col.axis=1, col.lab=1)
mtext("Average temperature", side=4, col=1, line=2)
grid(NA,NULL, lty = 1, col = colgrid) 
mtext("D",side=3, adj=0, font=2)

# Fig 2 - save ------------------------------------------------------------
dev.copy(pdf, paste0("output/fig2_modelOutputs", run.name,".pdf"),  12, 6)
  dev.off()

# Figure 3 - introduction dynamics ----------------------------------------
if(output_simulations==T){
  thetaMax <- thetatab[sim_liktab==max(sim_liktab),] %>% summarise_all(~median(.)) %>% as.numeric()
  thetaMax <- setNames(thetaMax, names(thetatab))
  thetaInitMax <- theta_inittab[sim_liktab==max(sim_liktab), ] %>% 
    summarise_all(~median(.)) %>% 
    as.numeric()
  thetaInitMax <- setNames(thetaInitMax, names(theta_inittab))
  thetaInitMax <- c(thetaInitMax, "c_init"=0, "cd_init"=0)
  names(thetaMax)[names(thetaMax)=="Inf."] <- "Inf"

  zikv_simulation <- function(change_startdate){
    ## loop through simulation
    thetaMax[["intro_mid"]] <- change_startdate
    sim <- Deterministic_modelR_final_DENVimmmunity(
      theta=c(thetaMax, theta_denv), 
      thetaInitMax, 
      seroposdates=seroposdates, 
      include.count=T
      )
    return(list(zika_trace = sim$C_trace, denv_trace = sim$CD_trace))
  }
  
  intros <- seq(365,365*2,90)
  zika_model_sims <- sapply(intros, function(xx){
    zikv_simulation(xx)$zika_trace})
  infections_series <- as.data.frame(zika_model_sims)
  names(infections_series) <- modelSt+intros
  infections_series$time_vals <- row(infections_series)[,1]*7
  infections_series$date_vals <- infections_series$time_vals + startdate
  infections_series$denv_vals <- zikv_simulation(1)$denv_trace
  
  short_infections_series <- infections_series %>%
    filter(date_vals >= as.Date("2013-10-01") & date_vals <= as.Date("2016-06-01"))
  seasonal_wave <- sapply(short_infections_series$time_vals, function(xx){seasonal_f(xx, amp = thetaMax[["beta_v_amp"]], mid=thetaMax[["beta_v_mid"]])})
  
  sims <- dim(short_infections_series)[2]-3
  pdf(here::here("output/fig3_simulations.pdf"), width = 8, height = 6)
  par(mfcol=c(sims, 1), mar=c(2,4,1,4)+0.1)
  ylimmax <- short_infections_series %>% pivot_longer(cols=starts_with("2")) %>% summarise(max(value)) %>% pull()
  ylimits <- c(0, round(ylimmax*1.1,0))
  for(ii in 1:sims){
    plot_date_series <- short_infections_series$date_vals
    ##
    #plot(short_infections_series$date_vals, seasonal_wave, col="gray85", yaxt="n", ylab="", xaxt="n", xlab="", type='l', ylim=c(0,2), bty="n", lwd = 1.25)
    #par(new=T)
    plot(short_infections_series$date_vals, short_infections_series[,ii], type="l", col=col1, ylab="", xlab="Date", ylim=ylimits, xaxt="n", bty="n", lwd = 2)
    mtext("Infections", side=2, line=2, col = col1)
    text(as.Date("2016-01-01"),ylimmax*(7/8), paste("Attack Rate","=",signif(sum(short_infections_series[,ii])/thetaMax[["npop"]],2)))
    lines(short_infections_series$date_vals, short_infections_series$denv_vals, lwd = 2 , col = col4)
    par(new=T)
    intros_series <- sapply(short_infections_series$time_vals, function(xx){intro_f(xx, mid = intros[ii], width = thetaMax[["intro_width"]], base = thetaMax[["intro_base"]])})
    plot(short_infections_series$date_vals, intros_series, type="l", yaxt="n", ylab="", xaxt="n", xlab="", col = col2, ylim=c(0,4), lty=2, bty="n", lwd = 1.25)
    axis.Date(side=1, at=seq.Date(min(plot_date_series), max(plot_date_series)+180, by = "6 months"), "months", format = "%b %Y")
    axis(side=4)
    mtext("Introductions", side=4, line=2, col = col2)
    mtext(LETTERS[ii], side=3, adj=0, font=2)
    grid(NA,NULL, lty = 1, col = colgrid) 
  }
  dev.off()
}
  

# MCMC diagnostics plots --------------------------------------------------
if(output_diagnostics==T){
  ## trace plots 
  labelN=1
  iiH=1
  parameters_to_estimate <- c(names(thetatab)[diag(var(thetatab))!=0 & !is.na(diag(var(thetatab)))])
  parameters_names <- stringr::str_replace_all(parameters_to_estimate, 
                                  pattern = c("r0" = "R0",
                                              "beta_h" = "Transmission rate",
                                              "rep" = "Reporting rate",
                                              "beta_base" = "Effect of clean-up campaign",
                                              "rec0" = "Initial immune proportion",
                                              "chi" = "Cross-protection from DENV3 outbreak",
                                              "epsilon" = "False positive rate in assay",
                                              "rho" = "Waning detectable antibody titres",
                                              "intro_mid" = "Midpoint of ZIKV introduction",
                                              "intro_base" = "Peak of introductions wave",
                                              "intro_width" = "Width of introductions wave")
                                  )
  names(parameters_names) <- parameters_to_estimate
  
  height_plot <- length(parameters_to_estimate)
  par(mfrow=c(ceiling(height_plot/2),4), mar  = c(2,3,1,3))
  ## function to split thetatab into chains and plot together
  plot_trace <- function(param){
    chain_length <- length(thetatab$npop)/m.tot
    d <- thetatab[, param]
    dsplit <- split(d, ceiling(seq_along(d)/chain_length))
    ylims <- c(max(0, min(d)-(0.25*min(d))), max(d)*1.25)
  plot(dsplit$`1`, type='l', col=col7, ylab="Parameter value", xlab="iteration", ylim=ylims,
       main=parameters_names[param])
  lines(dsplit$`2`, type='l', col=col3)
  lines(dsplit$`3`, type='l', col=col4)
  lines(dsplit$`4`, type='l', col=col5)
  
  hist(d, col=col1, xlab=param, main=paste0("ESS = ", signif(effectiveSize(thetatab[[param]]), 3)))
  }
    ## plot r0
    r0_chain <- apply(thetatab, 1, function(x) calculate_r0(th_in = as.list(x), sus = 342000, b_vary = 1)$r0_out)
    #thetatab$r0 <- r0_chain
  sapply(parameters_to_estimate, function(xx){plot_trace(xx)})
  dev.copy(pdf, paste0("output/fig5_diagnostics", virus, ".pdf"), width = length(parameters_to_estimate)*2, height=length(parameters_to_estimate))
    dev.off()
  
  ## correlation plot
  pdf(paste0("output/fig6_correlation",virus,".pdf"), width = length(parameters_to_estimate)*2, height = length(parameters_to_estimate)*2)
  n_param <- length(parameters_to_estimate)
  par(mfrow=c(n_param, n_param), mar  = c(2,3,1,3))
  
  thetatab0 <- thetatab %>% data.frame()
  sample.p <- sample(length(thetatab0$beta_h),1000,replace=T)
  thinner.theta <- thetatab0[sample.p,]
  thinner.theta <- thinner.theta %>% 
    select(all_of(parameters_to_estimate))
  
  for(ii in 1:n_param){
    for(jj in 1:n_param){
      if(ii<=jj){
        if(ii == jj){
          hist(thetatab0[[parameters_to_estimate[ii]]], main = parameters_names[parameters_to_estimate[ii]], xlab = NULL) 
        }else{
          plot(thinner.theta[[parameters_to_estimate[ii]]],thinner.theta[[parameters_to_estimate[jj]]],pch=19,cex=0.2, xlab="", ylab="",col="white",xaxt="n",yaxt="n",axes=F)
        }
      }else{
        plot(thinner.theta[[parameters_to_estimate[ii]]],thinner.theta[[parameters_to_estimate[jj]]],pch=19,cex=0.1, xlab="", ylab="")#, parameters_names[ii], ylab=parameters_names[jj])
        points(median(thinner.theta[[parameters_to_estimate[ii]]]),median(thinner.theta[[parameters_to_estimate[jj]]]),col="orange",cex=1.5,lwd=4)
      }
    }
  }
  dev.off()
}
  
# Parameter estimates -----------------------------------------------------
## parameter values table
c.text <- function(x,sigF=3){
  bp1=signif(c(median(x),quantile(x,0.025),quantile(x,0.975)),sigF)
  paste(bp1[1]," (",bp1[2],"-",bp1[3],")",sep="")
}

c.text.date <- function(x,sigF=3){
  bp1=signif(c(median(x),quantile(x,0.025),quantile(x,0.975)),sigF)
  paste(startdate+bp1[1]," (",startdate+bp1[2]," - ",startdate+bp1[3],")",sep="")
}

iiH=1
#load_posterior_1 <- load.posteriors(load.run.name=run.name, file.path="posterior_outputZ", iiH, mcmc.burn=mcmc.burn)
#  list2env(load_posterior_1,globalenv())

is.it.missing <- function(x){if(is.na(sum(x))){rep(0,max(picks))}else{x}}
beta_h1 <- is.it.missing(thetatab$beta_h )
Nsize <-   is.it.missing(thetatab$npop)
alpha_h <- is.it.missing(thetatab$Exp)
gamma <-   is.it.missing(thetatab$Inf.)
rep <-     is.it.missing(thetatab$rep )
repvol <- is.it.missing(thetatab$repvol )
amp <-     is.it.missing(thetatab$beta_v_amp)
mid <-     is.it.missing(thetatab$beta_v_mid)
inf0 <-  is.it.missing(thetatab$inf0)
rec0 <-  is.it.missing(thetatab$rec0)
omega <-  is.it.missing(thetatab$omega_d)
chi <-  is.it.missing(thetatab$chi)
epsilon <-  is.it.missing(thetatab$epsilon)
alpha <-  is.it.missing(thetatab$alpha)
rho <-  is.it.missing(thetatab$rho)
beta_base <-  is.it.missing(thetatab$beta_base)
beta_grad <-  is.it.missing(thetatab$beta_grad)
beta_mid <-  is.it.missing(thetatab$beta_mid)
zika_intro <-  is.it.missing(thetatab$intro_mid)
intro_width <-  is.it.missing(thetatab$intro_width)
intro_base <-  is.it.missing(thetatab$intro_base)

## total introductions
total_intro <- 4 * intro_width * intro_base # Integral of 4*base*exp(-(time-mid)/width)/(1+exp(-(time-mid)/width))^2 over -infty/infty is 4*base*width
median(total_intro)

# Calculate DIC
sim_likOut <- as.vector(sim_liktab)
pick.max = picks[sim_likOut[picks]==max(sim_likOut[picks], na.rm=T)][1]
loglik_theta_bar = sim_liktab[pick.max]
deviance.at.post.mean = -2*loglik_theta_bar 
effective.param = var(-2*sim_likOut)/2
dic.calc = deviance.at.post.mean + effective.param

param1 <- cbind(
  c.text(r0_post,3),
  c.text(rr_post,3),
  c.text(beta_h1[picks],2),
  c.text(1/alpha_h[picks],2),
  c.text(1/gamma[picks],2),
  c.text(rec0[picks],2),
  c.text(rep[picks],2),
  c.text(repvol[picks],2),
  c.text(amp[picks],2),
  c.text(mid[picks],2),
  c.text(omega[picks],2),
  c.text(chi[picks],2),
  c.text(epsilon[picks],2),
  c.text(alpha[picks],2),
  c.text(rho[picks],2),
  c.text(beta_base[picks],2),
  c.text(beta_grad[picks],2),
  c.text(beta_mid[picks],2),
  c.text.date(zika_intro[picks], 2),
  c.text(4 * intro_width[picks] * intro_base[picks], 4),
  c.text(intro_width[picks], 4),
  c.text(intro_base[picks], 4),
  max(sim_liktab),
  dic.calc
)

rownames(param1)=c(locationtab[iiH])
colnames(param1)=c(
  "R0", "RR",
  "beta_h",
  "Intrinsic incubation period",
  "Infectious period",
  "Initial immune", 
  "Reporting rate","Reporting dispersion",
  "Seasonal amplitude","Seasonal midpoint",
  "Cross immunity period","Cross protection",
  "False positives", 
  "Sensitivity", 
  "Waning Zika immunity",
  "base","grad",'mid',
  "Zika introduction date (mid)", "Number introduced", "Width introduced","Peak introduced",
  "max Likelihood",
  "DIC")
t(param1)
write.csv(t(param1),here::here("output",paste0("paramA",run.name,".csv")))
          