###===================================================================###
###                   Bayesian Project ZIP:                           ###
###       How habitat influences skunk detection and occurenc         ###
###===================================================================### 



#########################
#Set-Up
rm(list=ls()) #clears any objects from the environment
library(rjags)

setwd("...")
skunk.det <- as.matrix(read.csv("SR_Skunk_dets.csv")) #skunk dets from 30 sites and 185 occasions
shrub <-c(as.matrix(read.csv("SR_site_covs.csv"))) #whether each site has shrub (1) or grass (0) add c on the front so that it's a string, not a df

###############################
#specify data
n.site = nrow(skunk.det) #30
n.obs = ncol(skunk.det) #185 days

field.data <- list(n.site = n.site,  n.obs = n.obs, Y = skunk.det, shrub = shrub)

#write the model

model.string<- "
    model {

    # Priors
#with prior distribution, it goes mean, variance/precision for dnorm. Low precision allows a big spread. dnorm is the standard distribution to use for uniformed priors 

b0 ~ dnorm(0,.001)
b1 ~ dnorm(0,.001)
a0 ~ dnorm(0,.001) 
a1 ~ dnorm(0,.001)
    
    # Likelihood    

 #first we just write in each equation, then add the indexes, then make loops
  #make the loop for the site (j) first (the ecological process) then do the observation process (k) then put the equations into the loops. j/site equations go in j for loop, obs that have k also, go in the k loop
  
  for (j in 1:n.site){
z[j] ~ dbern(psi[j]) #bernoulli distribution because site is either occupied or not 
  logit(psi[j]) = b0+b1 * shrub[j] #shrub data doesn't change per visit, so only indexed by j 


  for(k in 1:n.obs){
Y[j,k] ~ dpois(lambda[j,k]) # use poisson because detections are a count bounded by time 
log(lambda[j,k]) = a0+a1 * shrub[j]



}
}

    # Derived quantities
#these are not key components in the model, but we could be interested in this

psi.shrub <- b0 + b1
psi.noshrub <- b0

lambda.shrub <- a0+a1
lambda.noshrub <- a0

	}"

skunk.model<-textConnection(model.string) #you can either save this as a text file on the computer then load it back in, this code keeps the model as text in the code


## Initial Values Creation
#we need to do this for the priors, this tells the MCMC which number to start guessing each of the priors and the N from, it's the initial value for each chain

inits <- function(){list( #says function here because the starting number is going to be drawn randomly from a distribution, it will draw a different number each time.
  b0 = rnorm(1,0,1), # 1= draw one time, 0=mean, 1=SD
  b1 = rnorm(1,0,1),
  a0 = rnorm(1,0,1), 
  a1 = rnorm(1,0,1),
  z = apply (skunk.det,1,max) +1)} #z can't be negative (bc it's occupancy), first value can't be 0 because it causes problems in the chain (although z can be 0 at some other point), so instead use the maximum number of values per site that you observed- you know there is at least that many. apply goes to counts, then goes to 1(rows) and takes the maximum value of each line then adds 1 to make sure there are no 0

inits()

##Select Parameters
# Specify which parameters we want to keep track of in the MCMC
params <- c("a0", "a1", "b0", "b1", "psi.shrub", "psi.noshrub", "lambda.shrub", "lambda.noshrub")


## MCMC settings
nc <- 3 #3 chains, people usually do 3
nb <- 1000 #number of burn-ins
ni <- 10000 #number of iterations
nt <- 1 # how often we will take a value from the chain. 1= take all the values, 2= discard half of the values, 10=discard 90% of values. don't thin unless it is essential to make output manageable


## Start Gibbs sampling (ie start the MCMC using the gibbs sampler (JAGS is only a gibbs sampler))
#get everything together, give it to jags

#this line builds the model, puts everything together
jags.mod <- jags.model(skunk.model,data=field.data ,n.chains=nc ,n.adapt =1000, inits=inits)
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

##Last step: backtransform the means so that we know our answer

#lambda.noshrub
exp(-4.77161) #0.008466738

#lambda.shrub
exp(-4.02417) #0.01787826

#psi.noshrub 
exp(-0.23865)/(exp(-0.23865)+1) #0.4406191

#psi.shrub      
exp(-0.06895)/(exp(-0.06895)+1) #0.4827693

