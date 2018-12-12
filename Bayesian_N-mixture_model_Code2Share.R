###===================================================================###
###              Case study: King Fire bird data                      ###
###                 Bayesian N-mixture model                          ###
###===================================================================### 

rm(list=ls())
#install.packages("rjags")
library(rjags)

setwd("...")

counts <- as.matrix(read.csv("HEWA_counts.csv"))
burn <- as.matrix(read.csv("Site_cov.csv"))
wind <- as.matrix(read.csv("Survey_cov.csv"))

# ===================== MODEL IN JAGS ======================================

## Specify data that will be used in the model

n.site = nrow(counts)
n.obs = ncol(counts)

field.data <- list(n.site = n.site,  n.obs = n.obs, Y = counts, burn = burn, wind = wind)

## Write model

model1.string <-"
    model {

    # Priors

	...
    
    # Likelihood    
    
	...

    # Derived quantities
	totalN <- ... # Estimate total population size across all sites

	}"

Nmixture_m<-textConnection(model1.string)

## Inits function

inits <- function(){list(
	...)}

## Parameters to estimate
params <- c(...)

## MCMC settings
nc <- ...
nb <- ...
ni <- ...
nt <- ... # dont thin unless it is essential to make output manageable

## Start Gibbs sampling
jags.mod <- jags.model(Nmixture_m ,data=... ,n.chains=... ,n.adapt =1000, inits=...)

chains <- coda.samples(jags.mod, params, ni, nt, n.burnin=nb)

length(chains)
head(chains[[1]])
nrow(chains[[1]])

## Results
summary(chains)

## Trace plots (see book p.173)
...

## Convrgence check via R-hat (see book p.173)
...

## Chains convergence graphic (see book p.171)
...

