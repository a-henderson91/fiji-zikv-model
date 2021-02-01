# - - - - - - - - - - - - - - - - - - - - - - - 
# Extract BEAST posterior and calculate parameters for informative prior on Intro of ZIKV
# Author: Alasdair Henderson
# github.com/a-henderson91/fiji-zikv-model
# - - - - - - - - - - - - - - - - - - - - - - - 
data_loc <- "beast/"
log_file <- "mono_2.log"
last_date <- as.Date("2017-02-20")

posterior <-
  read.delim(paste0(data_loc, log_file), comment.char = "#")
head(posterior)

nbIteration <- max(posterior$state)
burnIn <- 0.1 * nbIteration
burnInRow <- which(posterior$state >= burnIn)

fj.C.tmrca <- posterior$tmrca.fiji.C.[burnInRow]
summary(fj.C.tmrca)

## summary functions
c.text <- function(x, sigF = 3) {
  bp1 = signif(c(median(x), quantile(x, 0.025), quantile(x, 0.975)), sigF)
  paste(bp1[1], " (", bp1[2], "-", bp1[3], ")", sep = "")
}
c.text.date <- function(x, sigF = 3) {
  bp1 = signif(c(median(x), quantile(x, 0.025), quantile(x, 0.975)), sigF)
  paste(last_date - bp1[1],
        " (",
        last_date - bp1[3],
        " - ",
        last_date - bp1[2],
        ")",
        sep = "")
}

c.text(fj.C.tmrca, 3)
hist(fj.C.tmrca)
c.text.date(fj.C.tmrca * 365.25)
starting.times <- data.frame("x" = last_date - (fj.C.tmrca * 365.25))

# rescale the posterior to dates forward from 2013-10-11 ------------------
## currently dates back from 2017-02-20
adjustment <- as.numeric(last_date - as.Date("2013-01-01")) 
tmrca <- adjustment - (fj.C.tmrca * 365.25) 
hist(tmrca)

# fit empirical distribution ----------------------------------------------
fit <- MASS::fitdistr(tmrca, densfun="normal")  # we assume my_data ~ Normal(?,?)
hist(tmrca, pch=20, breaks=25, prob=TRUE, main="")
curve(dnorm(x, fit$estimate[1], fit$estimate[2]), col="red", lwd=2, add=T)

# Save mean and sd --------------------------------------------------------
intro_prior_mu <- fit$estimate[1] %>% as.numeric()
intro_prior_sigma <- fit$estimate[2] %>% as.numeric()

# plot MCMC convergence ---------------------------------------------------
ess_tmrca <- signif(coda::effectiveSize(posterior$tmrca.fiji.C.[burnInRow]),6)

png(here::here("output/BEAST_convergence.png"), res = 150, height = 6, width = 8, units = "in")
par(mfrow = c(2,1), mar = c(4,4,1,0.2))
plot(posterior$likelihood[burnInRow], type = "l", ylab = "Likelihood/prior dens.", col = col2, ylim = c(-4600, -4000))
lines(posterior$prior[burnInRow], col = col1)
legend("right", c("Prior", "Likelihood"), col=c(col1, col2), lwd=2, ncol = 1, cex = 1, pt.cex=1, inset=0.02, bty="n")
mtext("A", adj = 0, font = 2)

#plot(posterior$joint[burnInRow], type = "l", ylab = "Joint dist. (burnin)", col = col7a)
plot(posterior$tmrca.fiji.C.[burnInRow], type = "l", ylab = "tMRCA, Central Div.", col = col7)
mtext(paste0("ESS = ", ess_tmrca), adj = 1, font = 2)
mtext("B", adj = 0, font = 2)
dev.off()

tail(posterior)
