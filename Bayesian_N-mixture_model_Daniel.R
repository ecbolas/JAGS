###===================================================================###
###              Case study: King Fire bird data                      ###
###                 Bayesian N-mixture model                          ###
###===================================================================### 

rm(list=ls()) #clears any objects from the environment
#install.packages("rjags")
library(rjags) #Daniel used rjags because he thought it would be easier to extract values from the chain

setwd("...")

counts <- as.matrix(read.csv("HEWA_counts.csv")) #actual sample data of Hermit Warblers (54 sites were visited 3 times)
burn <-c(as.matrix(read.csv("Site_cov.csv"))) #whether each site was burned or not, abundance cov, these data are 1 or 0(unburned). add c on the front so that it's a string, not a df
wind <- as.matrix(read.csv("Survey_cov.csv")) #detection cov


# ===================== MODEL IN JAGS ======================================

## Specify data that will be used in the model

n.site = nrow(counts) #54
n.obs = ncol(counts) #3, number of visits

field.data <- list(n.site = n.site,  n.obs = n.obs, Y = counts, burn = burn, wind = wind) #combines all the data in one place, used burn and wind instead of random letters

## Write model

#write it in quotes below that becomes a text file to run in JAGS, which is similar to the language of R but not the same

#trying it first with weak priors= Low precision= big spread, so you want tau to be a small number
tau <- 1/1000 #.001, this is precision just putting it in manually for the priors

#when describing the model the order probably doesn't matter

model1.string <-"
    model {

    # Priors
b0 ~ dnorm(0,.001)
b1 ~ dnorm(0,.001)
a0 ~ dnorm(0,.001)
a1 ~ dnorm(0,.001)
   
  # Likelihood 
  #first we just write in each equation, then add the indexes, then make loops
  #make the loop for the site (j) first (the ecological process) then do the observation process (k) then put the equations into the loops. j/site equations go in j for loop, obs that have k also, go in the k loop
  
  for (j in 1:n.site){
  N[j] ~ dpois(lambda[j])
  log(lambda[j]) = b0+b1 * burn[j] #burn data doesn't change per visit, so only indexed by j
  
  for(k in 1:n.obs){
  Y[j,k] ~ dbin(p[j,k],N[j]) #refers to Y table, line j and columnn k
  logit(p[j,k]) = a0+a1 * wind[j,k]

  }
  }


    # Derived quantities
    #these are not key components in the model, but we could be interested in this

	totalN <-sum(N) # Estimate total population size across all sites

	}"

Nmixture_m<-textConnection(model1.string) #you can either save this as a text file on the computer then load it back in, this code keeps the model as text in the code

## Initial Values Creation
#we need to do this for the priors 
#this tells the MCMC which number to start guessing each of the priors and the N from, it's the initial value for each chain

inits <- function(){list( #says function here because the starting number is going to be drawn randomly from a distribution, it will draw a different number each time.
	b0 = rnorm(1,0,1), # 1= draw one time, 0=mean, 1=SD
	b1 = rnorm(1,0,1),
	a0 = rnorm(1,0,1),
	a1 = rnorm(1,0,1),
	N = apply (counts,1,max) +1)} #N can't be negative (bc it's abundance), first value can't be 0 because it causes problems in the chain (although N can be 0 at some other point), so instead use the maximum number of values per site that you observed- you know there is at least that many. apply goes to counts, then goes to 1(rows) and takes the maximum value of each line then adds 1 to make sure there are no 0

inits()

## Specify which parameters we want to keep track of in the MCMC
params <- c("a0", "a1", "b0", "b1", "N", "totalN")

## MCMC settings
nc <- 3 #3 chains, people usually do 3
nb <- 1000 #number of burn-ins
ni <- 10000 #number of iterations
nt <- 1 # how often we will take a value from the chain. 1= take all the values, 2= discard half of the values, 10=discard 90% of values. don't thin unless it is essential to make output manageable

## Start Gibbs sampling (ie start the MCMC using the gibbs sampler (JAGS is only a gibbs sampler))
#get everything together, give it to jags

#this line builds the model, puts everything together
jags.mod <- jags.model(Nmixture_m ,data=field.data ,n.chains=nc ,n.adapt =1000, inits=inits)
#plug in the model we specified, our data list, number of chains, n.adapt = ?, inits

#this is the actually doing of the MCMC to build the chains
chains <- coda.samples(jags.mod, params, ni, nt, n.burnin=nb) #coda.samples function you need the jags model object, other mcmc settings

length(chains)
head(chains[[1]]) #looking at beginning of first chain. Have one column for each paramter being estimated, each N, and totalN
nrow(chains[[1]]) #the number of iterations we specified

## Results
summary(chains)
#first tells you the stats for each parameter
#second part tells you the quantiles (desribes what the histogram looks like)

## Trace plots (see book p.173)
traceplot(chains[[1]]) #this looks terrible, probably because our priors started out so weak (b0 ~ dnorm(0,.001)) 

#alternative way to plot it
plot(1:ni,chains[[1]][,"a0"], type ="l" )




## Convergence check via R-hat (see book p.173)
...

## Chains convergence graphic (see book p.171)
...

