### Calculating 99% confidence intervals and 2nd derivatives from Generalized Additive Mixed Model fits
### for significantly nonlinear indicator-driver pairs.
### Generates output files of GAMM 2nd derivatives and ranges of bootstrapped 99% CI 
### for subsequent use in threshold identification
### Files can be used to find ranges where 2nd derivative CI is above or below 0 to identify threshold ranges;
### the best estimate of threshold point location is driver value where there is minimum or maximum in 2nd derivative along that range.

# Code by Shannon Hennessey - s.hennessey13@gmail.com
# Data from the Hawaii Monitoring and Reporting Collaborative (HIMARC)

rm(list = ls())

### Packages
library(reshape2)
library(data.table)
library(tidyverse)

### Set directories
setwd()
stats.dir <- paste0(getwd(),"/outputs/stats/")

### Set indicator, driver, and model type
resp.name <- "Grazer"
dri.name <- "rugosity_predicted"
mod <- "GAMM_binomial" # options are: "GAMM_binomial" or "GAMM_gamma_presence"

### Load files
dat <- read.csv(paste0(stats.dir, resp.name, "-", dri.name, "_", mod, "_trim.csv"), header = TRUE) # trimmed .csv file with indicator and driver data
gamm.data <- read.csv(paste0(stats.dir, resp.name, "-", dri.name, "_", mod, "_predicted.csv"), header = TRUE) # predicted model data frame
gamm.mat <- fread(paste0(stats.dir, resp.name, "-", dri.name, "_", mod, "_gamm.mat.csv"), header = TRUE) # bootstrap replicates

# make numeric sequence of indicator and driver values for models
myresponse <- dat[[resp.name]]
mydriver <- dat[[dri.name]]

# set CI
CIlevel <- 99
lowerQ <- (1 - (CIlevel/100))/2
upperQ <- 1 - lowerQ

# set to number of levels (integer) for each random effect
l.island <- length(unique(dat$island))  
l.Dataset <- length(unique(dat$Dataset))  
l.habitat <- length(unique(dat$habitat))  

# set parameters
sp.len <- 200 # Spline length
nb <- 1000  # Number of bootstrap replicates
dif1 <- 1 # First derivative
dif2 <- 2 # Second derivative

# average over the random effect combinations for the predicted gamm data
index <- rep(1:sp.len, l.island*l.Dataset*l.habitat)

dat.mean <- as.data.frame(aggregate(x = gamm.data[[1]], by = list(index), FUN = mean)$x)
dat.mean$driver <- gamm.data$driver[1:sp.len]
colnames(dat.mean) <- c('response','driver')


# for each bootstrap, average over the random effect combinations
mean.mat <- matrix(ncol = nb, nrow = sp.len)

for(i in 1:ncol(gamm.mat)){
  means <- aggregate(x = gamm.mat[[i]], by = list(index), FUN = mean)$x
  mean.mat[,i] <- means
}

rm(gamm.mat) # remove 'gamm.mat' to save memory


## Take the 1st and 2nd derivative of each bootstrap replicate, and calculate confidence intervals
# When 95% of the 1st or 2nd derivative replicates are greater or less than zero, 
# there is a significant trend (1st derivative), or threshold (2nd derivative).

## CI of GAMM bootstrap 
ci <- matrix(nrow = 2, ncol = sp.len) ## create a matrix to be filled with bootstrapped CI
rownames(ci) <- c("lower","upper")
for(j in 1:(sp.len)) {
  IC <- quantile(mean.mat[j,], c(lowerQ, upperQ), na.rm = TRUE)
  ci[,j] <- rbind(IC[1], IC[2])
}


## 1st Derivative Estimation
dif1.line <- diff(dat.mean$response, differences = 1) # actual 1st derivative estimate from original smoother
deriv.matrix.1 <- matrix(nrow = sp.len - dif1, ncol = nb) # create a matrix to for the 1st derivative estimates
for(k in 1:nb) {
  derivi <- mean.mat[,k]
  deriv.response <- diff(derivi, differences = dif1)
  driver.len <- length(mydriver) - dif1 
  deriv.object <- cbind(deriv.response, driver.len)
  deriv.matrix.1[,k] <- deriv.object[,1]
}

# CI of 1st derivative 
dif1.len <- sp.len - dif1
ci1 <- matrix(nrow = 2, ncol = dif1.len) ## create a matrix to be filled with bootstrapped CI
rownames(ci1) <- c("lower", "upper")

for(l in 1:dif1.len) {
  IC <- quantile(deriv.matrix.1[l,], c(lowerQ, upperQ), na.rm = TRUE)
  ci1[,l] <- rbind(IC[1], IC[2])
}


## 2nd Derivative Estimation 
dif2.line <- diff(dat.mean$response, differences = 2) # actual 2nd derivative estimate from original smoother
deriv.matrix.2 <- matrix(nrow = sp.len - dif2, ncol = nb) # create a matrix to for the 2nd derivative estimates
for(m in 1:nb) {
  derivi <- mean.mat[,m]
  deriv.response <- diff(derivi, differences = dif2)
  driver.len <- length(mydriver) - dif2
  deriv.object <- cbind(deriv.response, driver.len)
  deriv.matrix.2[,m] <- deriv.object[,1]
}

# CI of 2nd derivative 
dif2.len <- sp.len - dif2
ci2 <- matrix(nrow = 2, ncol = dif2.len) ## create a matrix to be filled with bootstrapped CI
rownames(ci2) <- c("lower","upper")

for(n in 1:dif2.len) {
  IC <- quantile(deriv.matrix.2[n,], c(lowerQ, upperQ), na.rm = TRUE)
  ci2[,n] <- rbind(IC[1], IC[2])
}


# make dataframe of CI 
ci.df <- as.data.frame(t(ci))
ci.df$driver <- dat.mean$driver
colnames(ci.df) <- c("lower", "upper", "driver")
ci.df <- ci.df[,c("driver","lower", "upper")]

# make dataframe of 1st deriv CI 
ci1.df <- as.data.frame(t(ci1))
ci1.df$derivative <- dif1.line 
ci1.df$driver <- dat.mean$driver[c(1:199)]
colnames(ci1.df) <- c("lower_ci1", "upper_ci1", "derivative_1", "driver")

# make dataframe of 2nd deriv CI to pinpoint where threshold is for driver
ci2.df <- as.data.frame(t(ci2))
ci2.df$derivative <- dif2.line 
ci2.df$driver <- dat.mean$driver[c(1:198)]
colnames(ci2.df) <- c("lower_ci2", "upper_ci2", "derivative_2", "driver")

ci.df <- left_join(dat.mean, ci.df, by = "driver")
ci.df <- left_join(ci.df, ci1.df, by = "driver")
ci.df <- left_join(ci.df, ci2.df, by = "driver")


# write file to use for threshold identification and plotting
write.csv(ci.df, paste0(stats.dir, resp.name, "-", dri.name, "_", mod, "_derivCI.csv"), row.names = FALSE)
