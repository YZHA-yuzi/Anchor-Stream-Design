##############################################################################
#                                                                            #
#                                                                            #
# NUMERICAL EXAMPLE: implement proposed estimators for a simulated data set  #
#                                                                            #
#                                                                            #
##############################################################################


## ---------------------------------------------------------------- ##
## NOTES: Please make sure the working directory is the folder 
## where the R script "FUN_AnchorStream.R" 
## and the simulated data set ("toydata.rda") are located at ##
## ---------------------------------------------------------------- ##

## load in packages
library(MCMCpack) 

## load in self-defined functions and read in simulated dataset 
source("FUN_AnchorStream.R")
load("toydata.rda") # dat.obs


## the simulated individual-level data 
head(dat.obs)

Ntot = 500; p2 = 0.1

## case counts estimation using individual-level data
re.counts <- AnchorStream_CaseCount(dat.obs = dat.obs,
                                    Ntot = Ntot, p2 = p2,
                                    num.post = 10000,
                                    seed = 1234,
                                    data.type = "individual",
                                    cellcounts.vec = NULL)

## case counts estimation using cell counts
### get cell counts as presented in Table 1 from individual-level data
nobs.bec <- get_nxbar_obs(dat.obs = dat.obs, Ntot = Ntot)$nvec
re.counts.1 <- AnchorStream_CaseCount(dat.obs = NULL,
                                      Ntot = Ntot, p2 = p2,
                                      num.post = 10000,
                                      seed = 1234,
                                      data.type = "aggregated",
                                      cellcounts.vec = nobs.vec)


## estimate continuous mean using individual-level data
re.continuous <- AnchorStream_Continuous(dat.obs = dat.obs,
                                         Ntot = Ntot, seed = 12345,
                                         nboot = 1000)






