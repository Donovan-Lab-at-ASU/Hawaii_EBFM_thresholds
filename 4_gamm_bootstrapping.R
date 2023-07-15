### Bootstrapping Generalized Additive Mixed Models for significantly nonlinear indicator-driver pairs
### Generates output file of bootstrapped GAMM fits used in CI estimation
# Code by Shannon Hennessey - s.hennessey13@gmail.com
# Data from the Hawaii Monitoring and Reporting Collaborative (HIMARC)

rm(list = ls())

### Packages
library(mgcv)

### Set directories
setwd()
stats.dir <- paste0(getwd(),"/outputs/stats/")

### Set indicator, driver, and model type
resp.name <- "Grazer"
dri.name <- "rugosity_predicted"
mod <- "GAMM_binomial" # options are: "GAMM_binomial" or "GAMM_gamma_presence"

### Load files
dat <- read.csv(paste0(stats.dir, resp.name, "-", dri.name, "_", mod, "_trim.csv"), header = TRUE) # trimmed .csv file with specific indicator and driver data


### Begin bootstrapping protocol: 
## 1) take a random sample of the data 1000 times
## 2) fit GAMM to each bootstrap replicate

## Set up GAMMs
sp.len <- 200 # Spline length
samp.n <- 3000 # number of samples to take for bootstrap
nb <- 1000  # Number of bootstrap replicates
ks <- 4 # constrain model 'knots' to 4 to prevent overfitting
set.seed(123)

# set to number of levels (integer) for each random effect -- used to generate new data for model prediction
l.island <- length(unique(dat$island))  
l.Dataset <- length(unique(dat$Dataset))  
l.habitat <- length(unique(dat$habitat))  

# Sample data 1000 times:
bootSample <- list()
for (i in seq_len(nb*2)){
  bootSample[[i]] <- dat[sample(seq_len(nrow(dat)), nrow(dat), replace=TRUE), c(dri.name, resp.name, "island", "Dataset", "habitat")]
}

# make matrix to be populated by each bootstrap replicate
gamm.mat <- matrix(ncol = nb, nrow = sp.len*l.island*l.Dataset*l.habitat)
j = 1

# fit GAMM for each bootstrap replicate
for(i in 1:(nb*2)) {
  dat.i <- as.data.frame(bootSample[[i]][1:5])
  dri.i <- dat.i[,1]
  resp.i <- dat.i[,2]
  gamm.i <- NULL
  
  if(mod == "GAMM_binomial"){
    try(gamm.i <- gam(resp.i ~ s(dri.i, bs= "tp", k = ks) +
                        s(island, bs = 're') +
                        s(Dataset, bs = 're') +
                        s(habitat, bs = 're'), 
                      optimMmethod="GCV.Cp", 
                      se = T, data = dat.i,
                      family=binomial(link = 'logit'), silent = TRUE))
  }
  
  if(mod == "GAMM_gamma_presence"){
    try(gamm.i <- gam(resp.i ~ s(dri.i, bs= "tp", k = ks) +
                        s(island, bs = 're') +
                        s(Dataset, bs = 're') +
                        s(habitat, bs = 're'),
                      optimMmethod="GCV.Cp", 
                      se = T, data = dat.i,
                      family=Gamma(link='log')), silent=TRUE)
  }
  
  if (length(gamm.i) != 0){
    newdata.i <- NULL
    try(newdata.i <- data.frame(dri.i = c(rep(rep(seq(from = min(dat[,dri.name]), 
                                                      to = max(dat[,dri.name]), 
                                                      length.out = sp.len)), (l.island*l.Dataset*l.habitat))), 
                                island = c(rep(rep(as.character(unique(dat$island)),
                                                   each = sp.len), 
                                               l.Dataset*l.habitat)),
                                Dataset = c(rep(rep(as.character(unique(dat$Dataset)),
                                                    each = sp.len*l.island), 
                                                l.habitat)),
                                habitat = c(rep(as.character(unique(dat$habitat)),
                                                each = sp.len*l.island*l.Dataset))),
        silent = TRUE)
    
    # if bootstrap doesn't have a certain level of random effect in sampled data, skip to next iteration
    if(length(newdata.i) == 0) {
      next
    }
    
    # append to gamm.mat
    try(gamm.mat[,j] <- predict.gam(gamm.i, newdata = newdata.i, type = "response", se.fit = T)$fit,
        silent = TRUE)
    
    # stop once you reach the desired number of bootstrap replicates
    if(j > nb){
      break
    }
    
    print(j)
    j = j+1
  }
}  # end bootstrap loop


## write gamm.mat to .csv, will be used to identify confidence intervals
write.csv(gamm.mat, paste0(stats.dir, resp.name, "-", dri.name, "_", mod, "_gamm.mat.csv"), row.names = FALSE)
