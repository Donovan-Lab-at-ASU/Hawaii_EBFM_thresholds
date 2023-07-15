### Predicting Generalized Additive Mixed Models for significantly nonlinear indicator-driver pairs
### Generates output file of predicted model data that is used in CI estimation
# Code by Shannon Hennessey - s.hennessey13@gmail.com
# Data from the Hawaii Monitoring and Reporting Collaborative (HIMARC)

rm(list = ls())

### Packages
library(mgcv)

### Set directories
setwd()
stats.dir <- paste0(getwd(),"/outputs/stats/")

### Load data
dat <- read.csv("https://raw.githubusercontent.com/Donovan-Lab-at-ASU/Hawaii_herbivore_management/main/data/donovan_etal_HIMARCdata_for_herbivore_publication.csv") # read in .csv file with columns of fish taxonomic group biomass and drivers
dat[,c("island","Dataset","habitat")] <- lapply(dat[,c("island","Dataset","habitat")], as.factor) # set random effects/categorical variables as factors

# set indicator and driver
resp.name <- "Grazer"
dri.name <- "rugosity_predicted"

# subset data by indicator, driver, and any random effect columns of interest
dat.cut <- dat[colnames(dat) %in% c(resp.name, dri.name, "island", "Dataset", "habitat")]
dat.cut <- dat.cut[apply(dat.cut, 1, function(x)!any(is.na(x))),] # remove columns with NA values

### Set up GAMMs
sp.len <- 200 # Spline length
ks <- 4 # # constrain model 'knots' to 4 to prevent overfitting

# set to number of levels (integer) for each random effect -- used to generate new data for model prediction
l.island <- length(unique(dat$island))  
l.Dataset <- length(unique(dat$Dataset))  
l.habitat <- length(unique(dat$habitat))  

# set model type (based on prior nonlinearity screening criteria)
mod <- "GAMM_binomial" # options are: "GAMM_binomial" or "GAMM_gamma_presence"

## re-fit GAMM to generate predicted data to use in CI estimation
# make column of transformed indicator presence/absence for fitting binomial model
dat.cut$binom <- ifelse(dat.cut[,resp.name] > 0, 1, 0) 

if(mod == "GAMM_binomial"){
  dat.binom <- dat.cut[,c("binom", dri.name, "island", "Dataset", "habitat")]
  colnames(dat.binom) <- c(resp.name, dri.name, "island", "Dataset", "habitat")
  
  # make numeric sequence of indicator and driver values for models
  myresponse <- dat.binom[[resp.name]]
  mydriver <- dat.binom[[dri.name]]
  
  try(gamm1 <- gam(myresponse ~ s(mydriver, bs= "tp", k = ks) + 
                     s(island, bs = 're') +
                     s(Dataset, bs = 're') +
                     s(habitat, bs = 're'),
                   optimMmethod = "GCV.Cp", 
                   se = T, data = dat.binom,
                   family = binomial(link = 'logit')))
}


if(mod == "GAMM_gamma_presence"){
  
  # make new data frame with presence-only data
  dat.presence <- dat.cut[which(dat.cut$binom == 1),]
  
  # re-define numeric sequence of indicator and driver values
  myresponse <- dat.presence[[resp.name]]
  mydriver <- dat.presence[[dri.name]]
  
  try(gamm1 <- gam(myresponse ~ s(mydriver, bs = "tp", k = ks) + 
                     s(island, bs = 're') +
                     s(Dataset, bs = 're') +
                     s(habitat, bs = 're'),
                   optimMmethod = "GCV.Cp", 
                   se = T, data = dat.presence,
                   family = Gamma(link = 'log')))
}


## generate new data to predict each model over range of driver values
new_data <- data.frame(mydriver = c(rep(rep(seq(from = min(mydriver), 
                                                to = max(mydriver), 
                                                length.out = sp.len)), (l.island*l.Dataset*l.habitat))),
                       island = c(rep(rep(as.character(unique(dat$island)),
                                          each = sp.len), 
                                      l.Dataset*l.habitat)),
                       Dataset = c(rep(rep(as.character(unique(dat$Dataset)),
                                           each = sp.len*l.island), 
                                       l.habitat)),
                       habitat = c(rep(as.character(unique(dat$habitat)),
                                       each = sp.len*l.island*l.Dataset)))

pred <- predict.gam(gamm1, newdata = new_data, type = "response", se.fit = T)

## bind predicted data
gamm.data <- data.frame(pred$fit, new_data)
colnames(gamm.data) <- c("response", "driver", "island", "Dataset", "habitat")

## save output files - GAMM model, predicted data, and trimmed dataframe for each model
saveRDS(gamm1, paste0(stats.dir, resp.name, "-", dri.name, "_", mod, ".rds"))
write.csv(gamm.data, paste0(stats.dir, resp.name, "-", dri.name, "_", mod, "_predicted.csv"), row.names = FALSE)

if(mod == "GAMM_binomial"){
  write.csv(dat.binom, paste0(stats.dir, resp.name, "-", dri.name, "_", mod, "_trim.csv"), row.names = FALSE)
}

if(mod == "GAMM_gamma_presence"){
  write.csv(dat.presence, paste0(stats.dir, resp.name, "-", dri.name, "_", mod, "_trim.csv"), row.names = FALSE)
}
