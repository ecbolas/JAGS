
rm(list=ls())

library(jagsUI)

# =========================================================================
# Linear Mixed-Effects Model: VIPER MASS - LENGTH FOR DIFFERENT POPULATIONS
#===========================================================================
# Does mass correlate with length in different viper populations?


v <- read.csv("visper.csv")
attach(v)

n.groups <- 56 # Number of populations
n.sample <- 10 # Number of vipers in each population
n <- n.groups * n.sample 

#equation:
# yi = aji + Bji * xi + ei
# === y ~ norm(tau, sigma) these two equations are the same
# yi = mass of snake i
# xi = length of snake i
#j = population
#aj ~ Normal(mu, sigma2) random effects for intercepts that are population-specific
#Bj ~ Normal(mu, sigma2) random effects for slopes that are population-specific
#ei ~Normal(0, sigma2) residual random effects

# ===================== MODEL IN FREQUENTIST =============================
#install.packages("lme4")
library('lme4') #package for linear mixed effects models
lme.fit1 <- lmer(mass ~ length + (1 | pop), REML = TRUE) #population as random effect, use REML, restricted maximum likelihood
lme.fit1

# ===================== MODEL IN JAGS ====================================
#equation:
# yi = aji + Bji * xi + ei
# === y ~ norm(tau, sigma) these two equations are the same



# Specify data that will be used in the model (this can be specified either before or after writing the model, but if done before can help decide how to write the model)

win.data <- list(mass = as.numeric(mass), pop = as.numeric(pop), length = length,
                 ngroups = max(as.numeric(pop)), n = n)

# Write model

#slightly different language then regular R that will read this as a whole sting of characters, you can comment stuff out like normal

cat("
    model {
    
    #Priors:
    
    #RANDOM INTERCEPT PER POPULATION
    #these means that alpha is calculated for every group
    for (i in 1:ngroups){
    alpha[i] ~ dnorm(mu.int, tau.int) #intercept
    } #this says that alpha is a paramter that follows a normal distribution with mu.int is       mean and tau.int is the variance
    
    #these two are the priors for alpha
    mu.int ~ dnorm(0, 0.001) # Mean hyperparameter for random intercept
    tau.int <- 1 / (sigma.int * sigma.int) #variance for random intercept, can't use sigma2       (variance) directly, you have to calculate tau from sigma
    sigma.int ~ dunif(0, 100) # SD hyperparameter for random intercept
    
    #COMMON SLOPE FOR ALL POPULATIONS, we know this so it is a prior
    beta ~ dnorm(0, 0.001) # Common slope
    
    #PRIORS FOR THE DATA (y)
    tau <- 1 / ( sigma * sigma) # Residual precision. can't use sigma2 (variance) directly,       you have to calculate tau from sigma
    sigma ~ dunif(0, 100) # Residual standard deviation
    


    #LIKELIHOOD
 
    
    for (i in 1:n) {
    mass[i] ~ dnorm(mu[i], tau) # The random variable
    mu[i] <- alpha[pop[i]] + beta* length[i] # Expectation
    }
    }
    ",fill = TRUE, file = "lmm.txt")

# Initial values for priors 
inits <- function(){list(alpha = rnorm(n.groups, 0, 2), beta = rnorm(1, 1, 1),
                         mu.int = rnorm(1, 0, 1), sigma.int = rlnorm(1), sigma = rlnorm(1))}

# Parameters to estimate
parameters <- c("alpha", "beta", "mu.int", "sigma.int", "sigma")

# MCMC settings
ni <- 2000 # iterations, anything >1000 good
nb <- 500 #burn-in (how many to delete in beginning)
nt <- 2 #thinning, 1 = not deleting anything
nc <- 3 # of chains


# Run the model in JAGS

out <- jags(win.data, inits, parameters, "lmm.txt", n.chains = nc,
            n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE) #it will estimate 3 chains, parallel true means each is calculated at the same time in a different part of the computer. parallel= false means they calculate them one at a time

# Inspect results
#good to simulate data first so you have coefficients to compare it to. data is simulated in LMM_JAGS_datasim
print(out, dig = 3)

traceplot(out)
traceplot(out, param = c("mu.int", "sigma.int", "beta"))

print(out, 2)

# Plot
par(mfrow = c(1,2))
plot(density(out$sims.list$mu.int), xlab="alpha.mean", ylab="Frequency", frame = F) 
abline(v = 230, col = "red", lwd = 3)

# Compare with input values??????
intercept.mean ; slope.mean ; intercept.sd ; slope.sd ; sd(eps)