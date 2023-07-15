### Plotting fitted GAMMs, the bootstrapped confidence interval, and the threshold range(s) and point estimate(s).
# Code by Shannon Hennessey - s.hennessey13@gmail.com
# Data from the Hawaii Monitoring and Reporting Collaborative (HIMARC)

rm(list = ls())

### Set directories
setwd()
figure.dir <- paste0(getwd(),"/outputs/plots/")
stats.dir <- paste0(getwd(),"/outputs/stats/")
code.dir <- "~/github/project_folder/"

### Load functions 
source(paste0(code.dir,'functions.R'))

### Set indicator, driver, model type
resp.name <- "Grazer"
dri.name <- "rugosity_predicted"
mod <- "GAMM_binomial" # options are: "GAMM_binomial" or "GAMM_gamma_presence"

### Load data
dat <- read.csv(paste0(stats.dir, resp.name, "-", dri.name, "_", mod, "_trim.csv"), header = TRUE) # trimmed .csv file with raw data
ci.mat <- read.csv(paste0(stats.dir, resp.name, "-", dri.name, "_", mod, "_derivCI.csv"), header = TRUE) # bootstrapped CI file

### Load file with all thresholds, one line for each significant indicator-driver pair
## Needs the following columns, in addition to model identifiers:
# 'thresh.min' - the minimum driver value of the lowest threshold range
# 'thresh.max' - the maximum driver value of the highest threshold range
# 'thresh1' - the driver value of the first threshold 'best estimate'
# 'thresh2' - the driver value of the second threshold 'best estimate', if applicable
thresholds <- read.csv(paste0(stats.dir, 'final_threshold_pairs.csv'), header = TRUE) 


### Plot final GAMM threshold response figure
gamm.threshold.fig(ci = ci.mat,
                   figure.dir = paste0(figure.dir, resp.name, "_", mod, "_threshold.png"),
                   myresponse = dat[[resp.name]],
                   mydriver = dat[[dri.name]],
                   resp.name = "Grazers",
                   dri.name = "Rugosity",
                   thresh.min = thresholds$thresh.min,
                   thresh.max = thresholds$thresh.max,
                   thresh1 = round(thresholds$thresh1, digits = 4),
                   thresh2 = round(thresholds$thresh2, digits = 4))
