
rm(list=ls())

library(jagsUI)

# =========================================================================
#               LMM: VIPER MASS - LENGTH FOR DIFFERENT POPULATIONS
#===========================================================================
# Does mass correlate with lenght in different visper populations?

setwd("C:/Users/Ana/Documents/PhD/Second chapter/Data/Examples")
v <- read.csv("visper.csv")
attach(v)

n.groups <- 56 # Number of populations
n.sample <- 10 # Number of vipers in each pop
n <- n.groups * n.sample 

# ===================== MODEL IN FREQUENTIST =============================

library('lme4')
lme.fit1 <- lmer(mass ~ length + (1 | pop), REML = TRUE)
lme.fit1

# ===================== MODEL IN JAGS ====================================

# Specify data that will be used in the model

win.data <- list(mass = as.numeric(mass), pop = as.numeric(pop), length = length,
                 ngroups = max(as.numeric(pop)), n = n)

# Write model
setwd("C:/Users/Ana/Documents/PhD/Second chapter/Data/Examples")

cat("
    model {
    
    # Priors:
    
    #RANDOM INTERCEPT PER POPULATION
    for (i in 1:ngroups){
    alpha[i] ~ dnorm(mu.int, tau.int) # Random intercepts
    }
    mu.int ~ dnorm(0, 0.001) # Mean hyperparameter for random intercepts
    tau.int <- 1 / (sigma.int * sigma.int)
    sigma.int ~ dunif(0, 100) # SD hyperparameter for random intercepts
    
    #COMMON SLOPE FOR ALL POPULATIONS
    beta ~ dnorm(0, 0.001) # Common slope
    
    tau <- 1 / ( sigma * sigma) # Residual precision
    sigma ~ dunif(0, 100) # Residual standard deviation
    
    # Likelihood
    
    for (i in 1:n) {
    mass[i] ~ dnorm(mu[i], tau) # The random variable
    mu[i] <- alpha[pop[i]] + beta* length[i] # Expectation
    }
    }
    ",fill = TRUE, file = "lmm.txt")

# Inits function
inits <- function(){list(alpha = rnorm(n.groups, 0, 2), beta = rnorm(1, 1, 1),
                         mu.int = rnorm(1, 0, 1), sigma.int = rlnorm(1), sigma = rlnorm(1))}

# Parameters to estimate
parameters <- c("alpha", "beta", "mu.int", "sigma.int", "sigma")

# MCMC settings
ni <- 2000
nb <- 500
nt <- 2
nc <- 3

# Run the model in JAGS

out <- jags(win.data, inits, parameters, "lmm.txt", n.chains = nc,
            n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

# Inspect results
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