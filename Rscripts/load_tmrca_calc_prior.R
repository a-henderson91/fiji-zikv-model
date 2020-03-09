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
fit <- fitdistr(tmrca, densfun="normal")  # we assume my_data ~ Normal(?,?)
hist(tmrca, pch=20, breaks=25, prob=TRUE, main="")
curve(dnorm(x, fit$estimate[1], fit$estimate[2]), col="red", lwd=2, add=T)

# Save mean and sd --------------------------------------------------------
intro_prior_mu <- fit$estimate[1] %>% as.numeric()
intro_prior_sigma <- fit$estimate[2] %>% as.numeric()
