### Fitting Generalized Additive Mixed Models for indicator-driver pairs.
### Generates output file of significantly nonlinear GAMMs to be used for subsequent analyses
# Code by Shannon Hennessey - s.hennessey13@gmail.com
# Data from the Hawaii Monitoring and Reporting Collaborative (HIMARC)

rm(list = ls())

### Packages
library(mgcv)

### Set directories
setwd()
stats.dir <- paste0(getwd(),"/outputs/stats/")
code.dir <- paste0(getwd(),"/code/")

### Load functions 
source(paste0(code.dir,'functions.R'))

### Load data
dat <- read.csv("https://raw.githubusercontent.com/Donovan-Lab-at-ASU/Hawaii_herbivore_management/main/data/donovan_etal_HIMARCdata_for_herbivore_publication.csv") # read in .csv file with columns of fish taxonomic group biomass and drivers
dat[,c("island","Dataset","habitat")] <- lapply(dat[,c("island","Dataset","habitat")], as.factor) # set random effects/categorical variables as factors

# create empty dataframe to append model diagnostics to
summary.df <- data.frame()
  
# set indicator and driver
resp.name <- "Grazer"
dri.name <- "rugosity_predicted"

# subset data by indicator, driver, and any random effect columns of interest
dat.cut <- dat[colnames(dat) %in% c(resp.name, dri.name, "island", "Dataset", "habitat")]
dat.cut <- dat.cut[apply(dat.cut, 1, function(x)!any(is.na(x))),] # remove columns with NA values


### Fit the GAMMs
ks <- 4 # constrain model 'knots' to 4 to prevent overfitting

## Fit GAMM: binomial (presence/absence data)
# make column of transformed indicator presence/absence for fitting binomial model
dat.cut$binom <- ifelse(dat.cut[,resp.name] > 0, 1, 0)

# make numeric sequence of indicator and driver values for models
myresponse <- dat.cut[[resp.name]]
mydriver <- dat.cut[[dri.name]]

try(gamm1 <- gam(binom ~ s(mydriver, bs = "tp",k = ks) + 
                   s(island, bs = 're') +
                   s(Dataset, bs = 're') +
                   s(habitat, bs = 're'),
                 optimMmethod = "GCV.Cp", 
                 se = T, data = dat.cut,
                 family = binomial(link = 'logit')))

# extract relevant model diagnostics
summary.gamm <- GAMM.summary(model = gamm1,
                             resp.name = resp.name, 
                             dri.name = dri.name, 
                             out.name = 'GAMM_binomial')

# bind to master data frame
summary.df <- rbind(summary.df, summary.gamm)


## Fit GAMM: gamma with log link (presence-only biomass data)
# make new data frame with presence-only data
dat.presence <- dat.cut[which(dat.cut$binom == 1),]

# re-define numeric sequence of indicator and driver values
myresponse <- dat.presence[[resp.name]]
mydriver <- dat.presence[[dri.name]]

try(gamm2 <- gam(myresponse ~ s(mydriver, bs = "tp", k = ks) + 
                   s(island, bs = 're') +
                   s(Dataset, bs = 're') +
                   s(habitat, bs = 're'),
                 optimMmethod = "GCV.Cp", 
                 se = T, data = dat.presence,
                 family = Gamma(link = 'log')))

# extract relevant model diagnostics
summary.gamm <- GAMM.summary(model = gamm2,
                             resp.name = resp.name, 
                             dri.name = dri.name, 
                             out.name = 'GAMM_gamma_presence')

# bind to master data.frame                               
summary.df <- rbind(summary.df, summary.gamm)

# rename columns
colnames(summary.df) <- c("Model", "Response", "Driver", "GCV", "dev.expl", "edf.driver","p.driver")

# select indicator-driver pairs that meet model nonlinearity significance criteria
best.mods <- summary.df[!is.na(summary.df$edf.driver) &
                         summary.df$edf.driver >= 2 &
                         summary.df$dev.expl >= 0.2 &
                         summary.df$p.driver < 0.05,] 

# save data frames for subsequent use in additional scripts
write.csv(summary.df, paste0(stats.dir, "gamm_models_diagnostics.csv"), row.names=FALSE)
write.csv(best.mods, paste0(stats.dir, "gamm_models_best.csv"), row.names=FALSE)

