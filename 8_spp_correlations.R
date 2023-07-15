### Script to calculate Pearson's correlation coefficients between species, families, and functional groups
# Code by Shannon Hennessey - s.hennessey13@gmail.com
# Data from the Hawaii Monitoring and Reporting Collaborative (HIMARC)

rm(list = ls())

### Packages
library(rstatix) # for correlation plots

### Load data 
## read in .csv file with columns of fish taxonomic group biomass and drivers
dat <- read.csv("https://raw.githubusercontent.com/Donovan-Lab-at-ASU/Hawaii_herbivore_management/main/data/donovan_etal_HIMARCdata_for_herbivore_publication.csv")

## Family-level correlation analysis
families <- c() # set vector of column names for inclusion in family-species comparisons

fam.df <- dat[, families] # subset data frame
names(fam.df) <- c() # set vector of names for plot axes

cor(fam.df) # family correlations
cor_plot(cor(fam.df)) # plot 


## Functional group-level correlation analysis
func.grps <- c() # set vector of column names for inclusion in functional group-species comparisons

func.df <- func.df[,func.grps] # subset data frame
names(func.df) <- c() # set vector of names for plot axes

cor(func.df) # functional group correlations
cor_plot(cor(func.df)) # plot 

