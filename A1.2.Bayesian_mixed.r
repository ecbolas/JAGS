
# =========================================================================
#               LMM: VIPER MASS - LENGTH FOR DIFFERENT POPULATIONS (p.151)
#===========================================================================
# ===================== DATA SIMULATION ===============================

n.groups <- 56 # Number of populations
n.sample <- 10 # Number of vipers in each pop
n <- n.groups * n.sample # Total number of data points
pop <- gl(n = n.groups, k = n.sample) # Indicator for population

#We directly normalize covariate length to avoid trouble with WinBUGS.
# Body length (cm)
original.length <- runif(n, 45, 70)
mn <- mean(original.length)
sd <- sd(original.length)
cat("Mean and sd used to normalise.original length:", mn, sd, "\n\n")
length <- (original.length-mn) / sd 
hist(length, col = "grey")

# We build a design matrix without intercept.
Xmat <- model.matrix(~ pop * length-1-length)
print(Xmat[1:21,], dig = 2) # Print top 21 rows

intercept.mean <- 230 # mu alpha
intercept.sd <- 20 # sigma alpha
slope.mean <- 60 # mu beta
slope.sd <- 30 # sigma beta
intercept.effects <- rnorm(n = n.groups, mean = intercept.mean, sd = intercept.sd)
slope.effects <- rnorm(n = n.groups, mean = slope.mean, sd = slope.sd)
all.effects <- c(intercept.effects, slope.effects) # Put them all together

# We assemble the measurements yi as before.
lin.pred <- Xmat[,] %*% all.effects # Value of lin.predictor
eps <- rnorm(n = n, mean = 0, sd = 30) # residuals
mass <- lin.pred + eps # response lin.pred + residual
hist(mass, col = "grey") # Inspect what we've created

library("lattice")
xyplot(mass ~ length | pop)

# ===================== MODEL IN FREQUENTIST =============================

library('lme4')
lme.fit1 <- lmer(mass ~ length + (1 | pop), REML = TRUE)
lme.fit1

# ===================== MODEL IN JAGS ===============================

library(jagsUI)

# Save data for exercise in class
d <- cbind(pop,mass,length)
colnames(d)[2] <- "mass"
setwd("C:/Users/Ana/Documents/PhD/Second chapter/Data/Examples")
write.csv(d,"visper.csv")

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

# Compare with input values
intercept.mean ; slope.mean ; intercept.sd ; slope.sd ; sd(eps)

# Plot
par(mfrow = c(1,2))
plot(density(out$sims.list$mu.int), xlab="alpha.mean", ylab="Frequency", frame = F) 
abline(v = 230, col = "red", lwd = 3)


# ----
# =========================================================================
#              GLMM: RED-BACKED SHRIKE ABUNDANCE IN DIFFERENT SITES (p.203)
#===========================================================================
# Pair counts over 30 years available in each of the 16 populations
# ===================== DATA SIMULATION ==================================

#We generate data under the random-coefficients model without correlation
#between the intercepts and slopes.
n.groups <- 16
n.years <- 30
n <- n.groups * n.years
pop <- gl(n = n.groups, k = n.years)
#We standardize the year covariate to a range from zero to one.
original.year <- rep(1:n.years, n.groups)
year <- (original.year-1)/29
#We build a design matrix without the intercept.
Xmat <- model.matrix(~ pop * year - 1 - year)
print(Xmat[1:91,], dig = 2) # Print top 91 rows
#Next, we draw the intercept and slope parameter values from their normal
#distributions, whose hyperparameters we need to select.
intercept.mean <- 3 # Choose values for the hyperparams
intercept.sd <- 1
slope.mean <- -2
slope.sd <- 0.6
intercept.effects <- rnorm(n = n.groups, mean = intercept.mean, sd = intercept.sd)
slope.effects <- rnorm(n = n.groups, mean = slope.mean, sd = slope.sd)
all.effects <- c(intercept.effects, slope.effects) # Put them all together

#We assemble the counts Ci by first computing the linear predictor, then
#exponentiating it and finally adding Poisson noise. Then, we look at the data.

lin.pred <- Xmat[,] %*% all.effects # Value of lin.predictor
C <- rpois(n = n, lambda = exp(lin.pred)) # Exponentiate and add Poisson noise
hist(C, col = "grey") # Inspect what we've created

#We use a lattice graph to plot the shrike counts against time for each population (Fig. 16.2).
library("lattice")
xyplot(C ~ original.year | pop, ylab = "Red-backed shrike counts", xlab = "Year")



# ===================== MODEL IN FREQUENTIST ===============================

# We assume that each population has a specific trend, i.e., that both slopes 
#and intercepts are independent random variables governed by common hyperparameters.

library('lme4')
glmm.fit <- glmer(C ~ year + (1 | pop) + ( 0 + year | pop), family = poisson)
summary(glmm.fit) # Inspect results


# ===================== MODEL IN JAGS ======================================

# Define model
library(jagsUI)

# Specify data that will be used in the model

win.data <- list(C = C, pop = as.numeric(pop), year = year, n.groups = n.groups, n = n)

# Write model
setwd("C:/Users/Ana/Documents/PhD/Second chapter/Data/Examples")

cat("
    model {

    # Priors
    for (i in 1:n.groups){
    alpha[i] ~ dnorm(mu.int, tau.int) # Intercepts
    beta[i] ~ dnorm(mu.beta, tau.beta) # Slopes
    }
    mu.int ~ dnorm(0, 0.001) # Hyperparameter for random intercepts
    tau.int <- 1 / (sigma.int * sigma.int)
    sigma.int ~ dunif(0, 10)
    mu.beta ~ dnorm(0, 0.001) # Hyperparameter for random slopes
    tau.beta <- 1 / (sigma.beta * sigma.beta)
    sigma.beta ~ dunif(0, 10)

    # Poisson likelihood
    for (i in 1:n) {
    C[i] ~ dpois(lambda[i])
    lambda[i] <- exp(alpha[pop[i]] + beta[pop[i]]* year[i])
    }
    }
    ", fill = TRUE, file = "glmm.txt")

# Eg: intercept.effects[pop[1]]

# Inits function
inits <- function(){ list(alpha = rnorm(n.groups, 0, 2), beta = rnorm(n.groups, 0, 2), 
                          mu.int = rnorm(1, 0, 1), sigma.int = rlnorm(1), mu.beta = rnorm(1, 0, 1),
                          sigma.beta = rlnorm(1))}
# Parameters to estimate
parameters <- c("alpha", "beta", "mu.int", "sigma.int", "mu.beta", "sigma.beta")
# MCMC settings
ni <- 2000
nb <- 500
nt <- 2
nc <- 3
# Start Gibbs sampling
out <- jags(win.data, inits, parameters, "glmm.txt", n.chains = nc,
            n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

print(out, dig = 2)

# ----
# =========================================================================
#   BINOMIAL MIXTURE GLMM: LIZARDS ABUNDANCE FROM DIFFERENT SURVEYS AND SITES
#===========================================================================
# ===================== DATA SIMULATION ====================================

# 1. Description of the biological process 

# We have surveyed 200 sites. Vegetation (quadratic) is affecting abundance
n.site <- 200
vege <- sort(runif(n = n.site, min = -1.5, max = 1.5))

#We construct the relationship between vegetation density and abundance

alpha.lam <- 2 # Intercept
beta1.lam <- 2 # Linear effect of vegetation
beta2.lam <- -2 # Quadratic effect of vegetation
lam <- exp(alpha.lam + beta1.lam * vege + beta2.lam * (vege^2))
par(mfrow = c(2,1))
plot(vege, lam, main ="", xlab = "", ylab = "Expected abundance", las = 1)
N <- rpois(n = n.site, lambda = lam)
table(N) # Distribution of abundances across sites
sum(N > 0) / n.site # Empirical occupancy
table(N)# Distribution of abundances across sites
plot(vege, N, main = "", xlab = "Vegetation cover", ylab = "Realized abundance")
points(vege, lam, type = "l", lwd = 2)

# 2. Description of the observation process

par(mfrow = c(2,1))
alpha.p <- 1 # Intercept
beta.p <- -4 # Linear effect of vegetation
det.prob <- exp(alpha.p + beta.p * vege) / (1 + exp(alpha.p + beta.p * vege))
plot(vege, det.prob, ylim = c(0,1), main = "", xlab = "", ylab = "Detection probability")

# Just for fun, see expected lizard count at each site in relation to vegetation cover
# The expected count (Mu) at a site i is given by True abundance (Ni) * p

expected.count <- N * det.prob
plot(vege, expected.count, main = "", xlab = "Vegetation cover", ylab = "Apparent
     abundance", ylim = c(0, max(N)), las = 1)
points(vege, lam, type = "l", col = "black", lwd = 2) # Truth

# 3. Simulate the data (What we would collect in the field)
# 3.1. Simulate the counts

R <- n.site
T <- 3 # Number of replicate counts at each site
y <- array(dim = c(R, T))
for(j in 1:T){
  y[,j] <- rbinom(n = n.site, size = N, prob = det.prob) # Size is the number you "through the coin"?
}

# Just for fun: Derive occupancy from abundance (the sum of occupied sites with N>0)

sum(apply(y, 1, sum) > 0) # Apparent distribution (proportion occupied sites)
sum(N > 0) # True occupancy

# Stack the replicated counts on top of each other for a vertical data
#format (i.e., convert the matrix to a vector)
C <- c(y)

# 3.2. Simulate, for the observation model, a site index and a vegetation covariate that
# have the same length as C (p suffix in variable name)
site = 1:R # ‘Short’ version of site covariate
site.p <- rep(site, T) # ‘Long’ version of site covariate
vege.p <- rep(vege, T) # ‘Long’ version of vegetation covariate
cbind(C, site.p, vege.p) # Check that all went right

# For fun: I think this would be the result when taking only the max count to model the effect
# of vegetation (without taking into account p)
max.count <- apply(y, 1, max)
naive.analysis <- glm(max.count ~ vege + I(vege^2), family = poisson)
summary(naive.analysis)
lin.pred <- naive.analysis$coefficients[1] + naive.analysis$coefficients[2] *
  vege + naive.analysis$coefficients[3] * (vege*vege)
# Plot comparing results and naive analysis
par(mfrow = c(1,1))
plot(vege, max.count, main = "", xlab = "Vegetation cover", ylab = "Abundance or count", 
     ylim = c(0,max(N)), las = 1)
points(vege, lam, type = "l", col = "black", lwd = 2)
points(vege, exp(lin.pred), type = "l", col = "red", lwd = 2)

# ===================== MODEL IN JAGS ======================================
# Define model
library(jagsUI)

# Specify data that will be used in the model

R = dim(y)[1]
n = dim(y)[1] * dim(y)[2]
vege2 = (vege * vege)
win.data <- list(R = R, vege = vege, vege2 = vege2, n = n, C = C, site.p = site.p,
                vege.p = vege.p)

# Write model
setwd("C:/Users/Ana/Documents/PhD/Second chapter/Data/Examples")
cat("
    model {

    # Priors

    alpha.lam ~ dnorm(0, 0.1)
    beta1.lam ~ dnorm(0, 0.1)
    beta2.lam ~ dnorm(0, 0.1)
    alpha.p ~ dnorm(0, 0.1)
    beta.p ~ dnorm(0, 0.1)

    # Likelihood

    # Biological model for true abundance
    for (i in 1:R) { # Loop over R sites
    N[i] ~ dpois(lambda[i])
    log(lambda[i]) <- alpha.lam + beta1.lam * vege[i] + beta2.lam * vege2[i]
    }
    # Observation model for replicated counts
    for (i in 1:n) { # Loop over all n observations
    C[i] ~ dbin(p[i], N[site.p[i]])
    logit(p[i]) <- alpha.p + beta.p * vege.p[i]
    }

    # Derived quantities
    totalN <- sum(N[]) # Estimate total population size across all sites

    }",fill = TRUE, file = "BinMix_pois.txt")

# Inits function
#clever starting values for the latent states (the Ni’s) are essential. 
#We use the maximum count at each site as a first guess of what N might be and 
#add 1 to avoid zeros. WinBUGS cannot use zeroes for the N value for the binomial 
#and will crash if you initialize the model with zeroes.
Nst <- apply(y, 1, max) + 1
inits <- function(){list(N = Nst, alpha.lam = rnorm(1, 0, 1), beta1.lam = rnorm(1, 0, 1), 
                         beta2.lam = rnorm(1, 0, 1), alpha.p = rnorm(1, 0, 1), beta.p = rnorm(1, 0, 1))}
# Parameters to estimate
params <- c("N", "totalN", "alpha.lam", "beta1.lam", "beta2.lam", "alpha.p",
           "beta.p")
# MCMC settings
nc <- 3
nb <- 200
ni <- 1200
nt <- 5
# Start Gibbs sampling
out <- jags(win.data, inits, params, "BinMix_pois.txt", n.chains = nc,
            n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)


# Let’s now compare the known true values with what the analysis has recovered.
cat("\n *** Our estimate of truth *** \n\n")
print(out, dig = 2)
cat("\n *** Compare with known truth ***\n\n")
alpha.lam ; beta1.lam ; beta2.lam ; alpha.p ; beta.p
sum(N) # True total population size across all sites
sum(apply(y, 1, max)) # Sum of site max counts
> cat("\n *** Our estimate of truth *** \n\n")
