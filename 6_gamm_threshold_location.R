### Identifying threshold ranges and point estimates from the 2nd derivative of bootstrapped Generalized Additive Mixed Model 99% CIs.
### Script finds ranges where 2nd derivative CI is above or below 0 to identify threshold ranges;
### the best estimate of threshold point location is driver value where there is minimum or maximum in 2nd derivative along that range.
# Code by Shannon Hennessey - s.hennessey13@gmail.com
# Data from the Hawaii Monitoring and Reporting Collaborative (HIMARC)

rm(list = ls())

### Packages
library(reshape2)

### Set directories
setwd()
stats.dir <- paste0(getwd(),"/outputs/stats/")

### Set indicator, driver, and model type
resp.name <- "Grazer"
dri.name <- "rugosity_predicted"
mod <- "GAMM_binomial" # options are: "GAMM_binomial" or "GAMM_gamma_presence"

# load files 
dat <- read.csv("https://raw.githubusercontent.com/Donovan-Lab-at-ASU/Hawaii_herbivore_management/main/data/donovan_etal_HIMARCdata_for_herbivore_publication.csv") # read in .csv file with columns of fish taxonomic group biomass and drivers
ci.mat <- read.csv(paste0(stats.dir, resp.name, "-", dri.name, "_", mod, "_derivCI.csv"), header = TRUE) # bootstrapped CI files

# create empty vector to append identified thresholds to
thresholds <- c()

## determine the signs of the upper and lower CI along the driver range
# if TRUE TRUE, both upper and lower CI are positive, if FALSE FALSE, both are negative, NA if 0
ci.mat$signs <-  ifelse(ci.mat$lower_ci2 != 0 & ci.mat$upper_ci2 != 0, 
                        ci.mat$signs <- paste(ci.mat$lower_ci2 > 0, ci.mat$upper_ci2 > 0, sep = "_"),
                        ci.mat$signs <- NA) 

## set driver minimum, maximum, and range
driver.min <- min(ci.mat$driver)
driver.max <- max(ci.mat$driver)
driver.range <- driver.max - driver.min

## estimate maxima
# select rows where both the upper and lower CI are above zero, ie. local maximum
temp.max <- ci.mat[which(ci.mat$signs == "TRUE_TRUE"),]

# calculate pairwise distances in driver values between each point (if maximum exists).
# if there is a large pairwise distance for one row, that indicates there are multiple thresholds that are CI maxima
try(temp.max$diffs <- c(0,diff(temp.max$driver)), silent = TRUE)

## if there is a threshold that is a CI maximum, subset first (or only) CI maximum 
if(nrow(temp.max) > 0){
  max1 <- c()
  row.count.max <- c()
  
  # iterate through rows that have the same pairwise distance between rows
  for(i in 1:nrow(temp.max)){
    temp <- temp.max[i,]
    if(temp$diffs < mean(temp.max$diffs[-c(1)]) + 0.01){
      max1 <- rbind(max1, temp)
    } else{break}
  }
  row.count.max <- rbind(row.count.max, i) 
  
  ## extract information on criteria to evaluate threshold robustness
  max.dat1 <- c(min(max1$driver), # threshold range minimum
                max(max1$driver), # threshold range maximum
                max1$driver[which(max1$derivative_2 == max(max1$derivative_2, na.rm = TRUE))], # threshold point estimate
                (min(max1$driver)-driver.min)/driver.range, # % of driver range threshold minimum falls within
                (max(max1$driver)-driver.min)/driver.range, # % of driver range threshold maximum falls within
                (max1$driver[which(max1$derivative_2 == max(max1$derivative_2, na.rm = TRUE))] - driver.min)/driver.range, # % of driver range threshold point estimate falls within
                max(max1$lower_ci2)) # CI maximum difference from zero
  
  ## if there is a second threshold that is a CI maximum, repeat identification process
  # select the rest of the rows to subset the second threshold that is a CI maximum
  if(max(row.count.max) < nrow(temp.max)){
    max2 <- c()
    for(i in max(row.count.max):nrow(temp.max)){
      temp <- temp.max[i,]
      max2 <- rbind(max2, temp)
    }
    
    ## extract information on criteria to evaluate threshold robustness
    max.dat2 <- c(min(max2$driver), # threshold range minimum
                  max(max2$driver), # threshold range maximum
                  max2$driver[which(max2$derivative_2 == max(max2$derivative_2, na.rm = TRUE))], # threshold point estimate
                  (min(max2$driver)-driver.min)/driver.range, # % of driver range threshold minimum falls within
                  (max(max2$driver)-driver.min)/driver.range, # % of driver range threshold maximum falls within
                  (max2$driver[which(max2$derivative_2 == max(max2$derivative_2, na.rm = TRUE))] - driver.min)/driver.range, # % of driver range threshold point estimate falls within
                  max(max2$lower_ci2)) # 2nd derivative CI maximum difference from zero
  }
}


## minima
# select rows where both the upper and lower CI are below zero, ie. local minimum
temp.min <- ci.mat[which(ci.mat$signs == "FALSE_FALSE"),]

# calculate pairwise distances in driver values between each point (if minimum exists).
# if there is a large pairwise distance for one row, that indicates there are multiple thresholds that are CI minima
try(temp.min$diffs <- c(0,diff(temp.min$driver)), silent = TRUE)

## if there is a threshold that is a CI minimum, subset first (or only) CI minimum
if(nrow(temp.min) > 0){
  min1 <- c()
  row.count.min <- c()
  for(i in 1:nrow(temp.min)){
    temp <- temp.min[i,]
    if(temp$diffs < (mean(temp.min$diffs[-c(1)]) + 0.01)){
      min1 <- rbind(min1, temp)
    } else{break}
    row.count.min <- rbind(row.count.min, i)
  }
  
  ## extract information on criteria to evaluate threshold robustness
  min.dat1 <- c(min(min1$driver), # threshold range minimum
                max(min1$driver), # threshold range maximum
                min1$driver[which(min1$derivative_2 == min(min1$derivative_2, na.rm = TRUE))], # threshold point estimate
                (min(min1$driver)-driver.min)/driver.range, # % of driver range threshold minimum falls within
                (max(min1$driver)-driver.min)/driver.range, # % of driver range threshold maximum falls within
                (min1$driver[which(min1$derivative_2 == min(min1$derivative_2, na.rm = TRUE))] - driver.min)/driver.range, # % of driver range threshold point estimate falls within
                min(min1$upper_ci2)) # 2nd derivative CI minimum difference from zero
  
  ## if there is a second threshold that is a CI maximum, repeat identification process
  # select the rest of the rows to subset the second threshold that is a CI maximum
  if(max(row.count.min) < nrow(temp.min)){
    min2 <- c()
    for(i in (max(row.count.min)+1):nrow(temp.min)){
      temp <- temp.min[i,]
      min2 <- rbind(min2, temp)
    }    
    
    ## extract information on criteria to evaluate threshold robustness
    min.dat2 <-  c(min(min2$driver), # threshold range minimum
                   max(min2$driver), # threshold range maximum
                   min2$driver[which(min2$derivative_2 == min(min2$derivative_2, na.rm = TRUE))], # threshold point estimate
                   (min(min2$driver)-driver.min)/driver.range, # % of driver range threshold minimum falls within
                   (max(min2$driver)-driver.min)/driver.range, # % of driver range threshold maximum falls within
                   (min2$driver[which(min2$derivative_2 == min(min2$derivative_2, na.rm = TRUE))] - driver.min)/driver.range, # % of driver range threshold point estimate falls within
                   min(min2$upper_ci2)) # 2nd derivative CI minimum difference from zero
  }
}

## bind everything together for that indicator-driver model
thresh.temp <- c()

if(exists("max.dat1")){
  thresh.temp <- rbind(thresh.temp,  max.dat1)
  rm(max.dat1, max1)
}

if(exists("min.dat1")){
  thresh.temp <- rbind(thresh.temp,  min.dat1)
  rm(min.dat1, min1)
}

if(exists("max.dat2")){
  thresh.temp <- rbind(thresh.temp,  max.dat2)
  rm(max.dat2, max2)
}

if(exists("min.dat2")){
  thresh.temp <- rbind(thresh.temp,  min.dat2)
  rm(min.dat2, min2)
}

## add column names and additional overarching model identifiers
if(!is.null(thresh.temp)){
  thresh.temp <- as.data.frame(thresh.temp)
  names(thresh.temp) <- c("lower", "upper", "point_est", "lower%", "upper%", "point_est%", "CI_diff0")
  
  thresh.temp$response <- c(resp.name)
  thresh.temp$driver <- c(dri.name)
  thresh.temp$model <- c(mod)
  
  thresholds <- rbind(thresholds, thresh.temp) # append to overall threshold df
}

## append count of how many data points within threshold range
n.points <- c()

for(i in 1:nrow(thresholds)){
  temp <- thresholds[i,]
  response <-  temp$response
  driver <- temp$driver
  
  points <- dat[, c(response,driver)]
  points$binom <-  ifelse(dat[,response] > 0, 1, 0)
  
  if(temp$model == "GAMM_binomial"){
    points <- points[which(points$binom == 1),]
    temp.n <- sum(points[,driver] > temp$lower & points[,driver] < temp$upper , na.rm=T)
  }
  
  if(temp$model == "GAMM_gamma_presence"){
    temp.n <- sum(points[,driver] > temp$lower & points[,driver] < temp$upper , na.rm=T)
  }
  
  n.points <- cbind(n.points, temp.n)
}

thresholds$n.points <- as.numeric(n.points)

# write file
write.csv(thresholds, paste0(stats.dir, "thresholds.csv"), row.names = FALSE)

