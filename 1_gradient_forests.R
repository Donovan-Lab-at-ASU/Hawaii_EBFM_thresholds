### Script for Gradient Forest analysis
# Code by Shannon Hennessey - s.hennessey13@gmail.com
# Data from the Hawaii Monitoring and Reporting Collaborative (HIMARC)

rm(list = ls())

### Packages
## May need to install packages with this code, package no longer available for more recent R versions
# install.packages("extendedForest", version = "1.6.1", repos="http://R-Forge.R-project.org")
# install.packages("gradientForest", version = "0.1-32", repos="http://R-Forge.R-project.org")
library(extendedForest)
library(gradientForest)

### Set directories
setwd()
stats.dir <- paste0(getwd(),"/outputs/stats/")

### Load data
dat <- read.csv("https://raw.githubusercontent.com/Donovan-Lab-at-ASU/Hawaii_herbivore_management/main/data/donovan_etal_HIMARCdata_for_herbivore_publication.csv") # read in .csv file with columns of fish taxonomic group biomass and drivers


### Gradient Forest analysis
## set variables
response <- c() # vector of indicator column names

# vector of driver column names
drivers <- c('OTP_CHL_ANOM_F', 'OTP_CHL_LTM', 'OTP_CHL_ANOM_M',   # Productivity
             'OTP_WAV_ANOM_F', 'OTP_WAV_ANOM_M',                  # Waves
             'OTP_PAR_LTM',                                       # Light
             'OTP_SST_LTM', 'OTP_SST_STD',                        # Temperature
             'rugosity_predicted', 'depth_predicted')             # Habitat

set.seed(123)
maxLevel <- log2(0.368*nrow(dat)/2) # calculate 'maxLevel' parameter per suggestion from Ellis et al. 2012

# run Gradient Forest model
gf <- gradientForest(dat,
                     predictor.vars = drivers, 
                     response.vars = response,
                     maxLevel = maxLevel,
                     ntree = 500, transform = NULL, compact = T,
                     nbin = 201, corr.threshold = 0.5)

# save indicator and driver R2 importances
write.csv(gf$imp.rsq, paste0(stats.dir,"gf_r2_importance.csv"), row.names = FALSE)

# select drivers with R2 importance above cutoff significance value
top_drivers <- names(importance(gf)[which(importance(gf) >=0.01)])


# plot R2 weighted importance of drivers
plot(gf, plot.type = "O")

# plot species cumulative importances
plot(gf, plot.type="C", imp.vars = top_drivers,
     show.overall = F, legend = T, leg.posn = "topleft", 
     leg.nspecies = 15, cex.lab = 0.7, cex.legend = 0.4, 
     cex.axis = 0.6, line.ylab = 0.9, par.args = list(mgp = c(1.5, 0.5, 0),
                                                      mar = c(2.5, 1, 0.1, 0.5),
                                                      omi = c(0, 0.3, 0, 0)))
