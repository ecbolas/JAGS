install.packages("pscl")
library(pscl)


zeroinfl(formula, data, offset,dist = "poisson",link = "log")

##formula has 2 components
## count component, modeling number of records if species is present
## and zero inflation component, modeling whether species is present (occupancy)

## for example y~dist|1 is a model where the count component is influenced by the
## covariate "dist", and the zero inflation component has no covariates

## for example y~dist|dist is a model where the count component 
## and the zero inflation component are influenced by covariate "dist"


###generate some fake data

J=100 #number of sites
p0<-0.2 #percent sites not occupied
dist<-rnorm(J,0,1) ##covariate "dist"
alpha<-0.05 #log-linear intercept, photographic rate
beta<-0.5 #positive effect of "dist" on photographic rate

##make effort - this would be the number of days the camera worked
eff<-sample(1:5, J, replace=T) ##effort varies between 1 and 5


##calculate expected photographic rate for all J cameras, based on dist and effort (=offset)
lam<-exp(alpha+beta*dist + log(eff))


##create indicator whether site is occupied or not
z<-rbinom(J,1, (1-p0))

##effective expected photographic rate
lam.eff<-lam*z

##generate records - this would be your number of records for a species at each camera
y<-rpois(J, lam.eff)

###put everything in a data frame for analysis
dat<-data.frame(y=y, dist=dist, eff=eff)

##run zero inflated models
##note that effort is used as an offset and needs to be provided on the log-scale
mod0<-zeroinfl(y~1|1, data=dat, dist = "poisson", offset=log(eff)) #null model
mod1<-zeroinfl(y~dist|1, data=dat, dist = "poisson", offset=log(eff)) #"dist" on count data
mod2<-zeroinfl(y~dist|dist, data=dat, dist = "poisson", offset=log(eff)) #"dist" on count data and zero inflation

##compare by AIC
AIC(mod0, mod1, mod2)


###looking at model output
summary(mod1)

#[...]

#Count model coefficients (poisson with log link):
#            Estimate Std. Error z value Pr(>|z|)    
#(Intercept) 0.08498    0.07439   1.142    0.253    
#dist         0.52117    0.06618   7.875 3.41e-15 ***

#[Intercept and covariate coefficient for the number of records]

#Zero-inflation model coefficients (binomial with logit link):
#            Estimate Std. Error z value Pr(>|z|)  
#(Intercept)  -1.455      0.322  -4.519 6.22e-06 ***

#[use plogis(-1.455) to get probability that site is NOT occupied and 
# 1-plogis(-1.455) to get probability that site is occupied]

#---
#Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 

#Number of iterations in BFGS optimization: 11 
#Log-likelihood: -185.6 on 3 Df











