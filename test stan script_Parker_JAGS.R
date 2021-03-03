library(pivmet)
library(bayesmix)
library(bayesplot)
library(tidyverse)
library(R2jags) 
library(mixtools)
library(coda)
library(lattice)


#Resource used to work out how to extract loglikelihood and calc WAIC and LOO.
#https://discourse.mc-stan.org/t/calculating-the-likelihood-for-loo-in-a-binomial-mixture-model/4787/5

options(mc.cores = 4)
rstan_options(auto_write = TRUE)

#https://cran.r-project.org/web/packages/pivmet/vignettes/Relabelling_in_Bayesian_mixtures_by_pivotal_units.html

## Laterlity data: Parker et al. (2020)

### Data analysis

# read in data from csv
#LIs <- read.csv('https://osf.io/x5jd2/download',stringsAsFactors = F)
#We will replace with read from OSF, but better stick with LIs to ensure final version is the one we use 
allsum <- read.csv('./data/allsum_out.csv',stringsAsFactors = F, na.strings = "NA")
allsum <- allsum[!is.na(allsum$origID), ]



# remove those without day 2 to make LIs file - we will just use this
# reformat NA
#identify missing cases for time 2
day2cols<-c('Day2_CF_acc_LI','Day2_DL_acc_LI','Day2.RDT_RT_LI','Day2_finger_LI')
allsum$day2.complete <- ifelse(is.na(allsum$day2.complete), 0, allsum$day2.complete)
LIs <- allsum[allsum$day2.complete == 1,]

# remove
LIs_c <- LIs[LIs$day2.complete == 1,]

# now remove them from LI
LIs <- LIs[LIs$lexOut == "keep",]
# table of counts for footed and sigh
short_dat <- dplyr::select(LIs, "subject", "handedness", "footedness", "miles", "porta")
long_dat <- tidyr::gather(short_dat, task, response, footedness:porta, factor_key=TRUE)

#Add a category of handedness that distinguishes extreme left from other left-handers
#In preregistration, extreme left defined as EHI laterality -90 or less.
LIs$handed3 <- LIs$handedness
w<-which(LIs$index_EHI<(-89))
LIs$handed3[w]<-'Extreme.left'

myLI<-LIs$Day1_DL_acc_LI[complete.cases(LIs$Day1_DL_acc_LI)]
myLI2<-LIs$Day1_CF_acc_LI[complete.cases(LIs$Day1_CF_acc_LI)]
#==================================================================##

plot.normal.components <- function(mixture,component.number,...) {
  curve(mixture$lambda[component.number] *
          dnorm(x,mean=mixture$mu[component.number],
                sd=mixture$sigma[component.number]), add=TRUE, ...)
}

dnormalmix <- function(x,mixture,log=FALSE) {
  lambda <- mixture$lambda
  k <- length(lambda)
  # Calculate share of likelihood for all data for one component
  like.component <- function(x,component) {
    lambda[component]*dnorm(x,mean=mixture$mu[component],
                            sd=mixture$sigma[component])
  }
  # Create array with likelihood shares from all components over all data
  likes <- sapply(1:k,like.component,x=x)
  # Add up contributions from components
  d <- rowSums(likes)
  if (log) {
    d <- log(d) }
  return(d) }
loglike.normalmix <- function(x,mixture) {
  loglike <- dnormalmix(x,mixture,log=TRUE)
  return(sum(loglike))
}

#==================================================================##
# FREQUENTIST APPROACH via EM algorithm.
#==================================================================##

# Basic EM algorithm script using mixtools (https://www.stat.cmu.edu/~cshalizi/402/lectures/20-mixture-examples/lecture-20.pdf)

LI_fit_DL<-normalmixEM(myLI)
plot(LI_fit_DL, density=TRUE, cex.axis=1.4, cex.lab=1.4, cex.main=1.8,main2="Dichotic Listening", xlab2="LI")

LI_fit_CF<-normalmixEM(myLI2)
plot(LI_fit_CF, density=TRUE, cex.axis=1.4, cex.lab=1.4, cex.main=1.8,main2="Chimeric faces", xlab2="LI")



LI.DL.boot <- boot.comp(myLI,max.comp=10,mix.type="normalmix", maxit=400,epsilon=1e-2)

LI.CF.boot <- boot.comp(myLI2,max.comp=10,mix.type="normalmix", maxit=400,epsilon=1e-2)

#Dichotic Listening#

LI.DL.k2 <- normalmixEM(myLI,k=2,maxit=100,epsilon=0.01)
LI.DL.k3 <- normalmixEM(myLI,k=3,maxit=100,epsilon=0.01)
LI.DL.k4 <- normalmixEM(myLI,k=4,maxit=100,epsilon=0.01)

plot(hist(myLI,breaks=101),col="grey",border="grey",freq=FALSE,
     xlab="LI",main="Dichotic Listening (k=2)")
lines(density(myLI),lty=2)
sapply(1:2,plot.normal.components,mixture=LI.DL.k2)


plot(hist(myLI,breaks=101),col="grey",border="grey",freq=FALSE,
     xlab="LI",main="Dichotic Listening (k=3)")
lines(density(myLI),lty=2)
sapply(1:3,plot.normal.components,mixture=LI.DL.k3)

plot(hist(myLI,breaks=101),col="grey",border="grey",freq=FALSE,
     xlab="LI",main="Dichotic Listening (k=4)")
lines(density(myLI),lty=2)
sapply(1:4,plot.normal.components,mixture=LI.DL.k4)

#Chimeric faces#

LI.CF.k2 <- normalmixEM(myLI2,k=2,maxit=100,epsilon=0.01)
LI.CF.k3 <- normalmixEM(myLI2,k=3,maxit=100,epsilon=0.01)
LI.CF.k4 <- normalmixEM(myLI2,k=4,maxit=100,epsilon=0.01)

plot(hist(myLI2,breaks=101),col="grey",border="grey",freq=FALSE,
     xlab="LI",main="Chimeric faces (k=2)")
lines(density(myLI2),lty=2)
sapply(1:2,plot.normal.components,mixture=LI.CF.k2)


plot(hist(myLI2,breaks=101),col="grey",border="grey",freq=FALSE,
     xlab="LI",main="Chimeric faces (k=3)")
lines(density(myLI2),lty=2)
sapply(1:3,plot.normal.components,mixture=LI.CF.k3)

plot(hist(myLI2,breaks=101),col="grey",border="grey",freq=FALSE,
     xlab="LI",main="Chimeric faces (k=4)")
lines(density(myLI2),lty=2)
sapply(1:4,plot.normal.components,mixture=LI.CF.k4)

####################################################################
# Cross validation to check number of components.

CVfun<-function(data)
{
  n <- length(data)
  data.points <- 1:n
  data.points <- sample(data.points) # Permute randomly
  train <- data.points[1:floor(n/2)] # First random half is training
  
  test <- data.points[-(1:floor(n/2))] # 2nd random half is testing
  candidate.component.numbers <- 2:10
  loglikes <- vector(length=1+length(candidate.component.numbers))
  # k=1 needs special handling
  mu<-mean(data[train]) # MLE of mean
  sigma <- sd(data[train])*sqrt((n-1)/n) # MLE of standard deviation
  loglikes[1] <- sum(dnorm(data[test],mu,sigma,log=TRUE))
  for (k in candidate.component.numbers) {
    mixture <- normalmixEM(data[train],k=k,maxit=400,epsilon=1e-2)
    loglikes[k] <- loglike.normalmix(data[test],mixture=mixture)
  }
  
  plot(x=1:10, y=loglikes,xlab="Number of mixture components",
       ylab="Log-likelihood on testing data")
}

CVfun(myLI)
CVfun(myLI2)

#Pretty inconclusive!
#==================================================================##
#==================================================================##
#
# BAYESIAN APPROACH!!
#
#==================================================================##
#
# JAGS model

cat("
model	{
  for (i in 1:N) {
    y[i] ~ dnorm(mu[S[i]],tau[S[i]]);
    S[i] ~ dcat(eta[]);
  }
  for (j in 1:k) {
    mu[j] ~ dnorm(b0,B0inv);
    tau[j] ~ dgamma(nu0Half,nu0S0Half);
    sigma[j] <- 1/tau[j]
  }
  S0 ~ dgamma(g0Half,g0G0Half);
  nu0S0Half <- nu0Half * S0;
  
  eta[1:k] ~ ddirch(e[1:k]);
}
",file="JAGS_mod1.txt")


##########
# Make a function for generating initial guesses

init.generator1 <- function(){ list(
  mu = runif(2, -100,100),
  tau = 1/runif(2, 1,10),
  eta = runif(2, 0,1))
}

init.generator2 <- function(){ list(
  mu = runif(3, -100,100),
  tau = 1/runif(3, 1,10),
  eta = runif(3, 0,1))
}

init.generator3 <- function(){ list(
  mu = runif(4, -100,100),
  tau = 1/runif(4, 1,10),
  eta = runif(4, 0,1))
}

#==================================================================##
# Dichotic listening
#==================================================================##
# 2 component Gaussian mixture model via Stan
#==================================================================##
#
#######
# Package the data for JAGS

data.package.DL2 <- list(N=length(myLI), y=myLI, k=2,
                      b0=0, B0inv=0.1,
                      nu0Half=1,
                      g0Half=1,g0G0Half=1, e=c(1,1))


params.to.monitor <- c("eta","mu","tau")

jags.fit.DL2 <- jags(data=data.package.DL2,inits=init.generator1,parameters.to.save=params.to.monitor,n.iter=10000,model.file="JAGS_mod1.txt",n.chains = 4,n.burnin = 2000,n.thin=5 )

jagsfitDL2.mcmc <- as.mcmc(jags.fit.DL2) 

summary(jagsfitDL2.mcmc)

plot(jagsfitDL2.mcmc[,"eta[1]"])
plot(jagsfitDL2.mcmc[,"eta[2]"])
plot(jagsfitDL2.mcmc[,"mu[1]"])
plot(jagsfitDL2.mcmc[,"mu[2]"])
plot(jagsfitDL2.mcmc[,"tau[1]"])
plot(jagsfitDL2.mcmc[,"tau[2]"])

densityplot(jagsfitDL2.mcmc)

DIC_DL2 <- jags.fit.DL2$BUGSoutput$DIC



#==================================================================##
k <- 2
nMC <- 10000

res2_DL <- piv_MCMC(y = myLI, k = k, nMC = nMC, 
                    software = "rjags")
rel2_DL <- piv_rel(res2_DL)

png(filename = "Parker_univ2_DL_trace_jags.png", width = 10, height = 4, units = "in",res=300)
piv_plot(y=myLI, res2_DL, rel2_DL, par = c("mean"), type="chains")
dev.off()

png(filename = "Parker_univ2_DL_hist_jags.png", width = 10, height = 4, units = "in",res=300)
piv_plot(y=myLI, res2_DL, rel2_DL, par = c("mean"), type="hist")
dev.off()

#==================================================================##
# 3 component Gaussian mixture model via Stan
#==================================================================##

#######
# Package the data for JAGS

data.package.DL3 <- list(N=length(myLI), y=myLI, k=3,
                         b0=0, B0inv=0.1,
                         nu0Half=1,
                         g0Half=1,g0G0Half=1, e=c(1,1,1))



params.to.monitor <- c("eta","mu","tau")

jags.fit.DL3 <- jags(data=data.package.DL3,inits=init.generator2,parameters.to.save=params.to.monitor,n.iter=10000,model.file="JAGS_mod1.txt",n.chains = 4,n.burnin = 2000,n.thin=5 )

jagsfitDL3.mcmc <- as.mcmc(jags.fit.DL3) 

summary(jagsfitDL3.mcmc)

plot(jagsfitDL3.mcmc[,"eta[1]"])
plot(jagsfitDL3.mcmc[,"eta[2]"])
plot(jagsfitDL3.mcmc[,"mu[1]"])
plot(jagsfitDL3.mcmc[,"mu[2]"])
plot(jagsfitDL3.mcmc[,"tau[1]"])
plot(jagsfitDL3.mcmc[,"tau[2]"])

densityplot(jagsfitDL3.mcmc)

DIC_DL3 <- jags.fit.DL3$BUGSoutput$DIC


#==================================================================##
k <- 3
nMC <- 10000

res3_DL <- piv_MCMC(y = myLI, k = k, nMC = nMC, 
                    software = "rjags")
rel3_DL <- piv_rel(res3_DL)

png(filename = "Parker_univ3_DL_trace_jags.png", width = 10, height = 4, units = "in",res=300)
piv_plot(y=myLI, res3_DL, rel3_DL, par = c("mean"), type="chains")
dev.off()

png(filename = "Parker_univ3_DL_hist_jags.png", width = 10, height = 4, units = "in",res=300)
piv_plot(y=myLI, res3_DL, rel3_DL, par = c("mean"), type="hist")
dev.off()

#==================================================================##
# 4 component Gaussian mixture model via Stan
#==================================================================##

#######
# Package the data for JAGS

data.package.DL4 <- list(N=length(myLI), y=myLI, k=4,
                         b0=0, B0inv=0.1,
                         nu0Half=1,
                         g0Half=1,g0G0Half=1, e=c(1,1,1,1))

params.to.monitor <- c("eta","mu","tau")

jags.fit.DL4 <- jags(data=data.package.DL4,inits=init.generator3,parameters.to.save=params.to.monitor,n.iter=10000,model.file="JAGS_mod1.txt",n.chains = 4,n.burnin = 2000,n.thin=5 )

jagsfitDL4.mcmc <- as.mcmc(jags.fit.DL4) 

summary(jagsfitDL4.mcmc)

plot(jagsfitDL4.mcmc[,"eta[1]"])
plot(jagsfitDL4.mcmc[,"eta[2]"])
plot(jagsfitDL4.mcmc[,"mu[1]"])
plot(jagsfitDL4.mcmc[,"mu[2]"])
plot(jagsfitDL4.mcmc[,"tau[1]"])
plot(jagsfitDL4.mcmc[,"tau[2]"])

densityplot(jagsfitDL4.mcmc)

DIC_DL4 <- jags.fit.DL4$BUGSoutput$DIC

#==================================================================##
k <- 4
nMC <- 20000

res4_DL <- piv_MCMC(y = myLI, k = k, nMC = nMC, 
                    software = "rjags")
rel4_DL <- piv_rel(res4_DL)

png(filename = "Parker_univ4_DL_trace_jags.png", width = 10, height = 4, units = "in",res=300)
piv_plot(y=myLI, res4_DL, rel4_DL, par = c("mean"), type="chains")
dev.off()

png(filename = "Parker_univ4_DL_hist_jags.png", width = 10, height = 4, units = "in",res=300)
piv_plot(y=myLI, res4_DL, rel4_DL, par = c("mean"), type="hist")
dev.off()

#==================================================================##
#==================================================================##
# Chimeric faces
#==================================================================##
# 2 component Gaussian mixture model via Stan
#==================================================================##

#######
# Package the data for JAGS

data.package.CF2 <- list(N=length(myLI), y=myLI, k=2,
                         b0=0, B0inv=0.1,
                         nu0Half=1,
                         g0Half=1,g0G0Half=1, e=c(1,1))

params.to.monitor <- c("eta","mu","tau")

jags.fit.CF2 <- jags(data=data.package.CF2,inits=init.generator1,parameters.to.save=params.to.monitor,n.iter=10000,model.file="JAGS_mod1.txt",n.chains = 4,n.burnin = 2000,n.thin=5 )

jagsfitCF2.mcmc <- as.mcmc(jags.fit.CF2) 

summary(jagsfitCF2.mcmc)

plot(jagsfitCF2.mcmc[,"eta[1]"])
plot(jagsfitCF2.mcmc[,"eta[2]"])
plot(jagsfitCF2.mcmc[,"mu[1]"])
plot(jagsfitCF2.mcmc[,"mu[2]"])
plot(jagsfitCF2.mcmc[,"tau[1]"])
plot(jagsfitCF2.mcmc[,"tau[2]"])

densityplot(jagsfitCF2.mcmc)

DIC_CF2 <- jags.fit.CF2$BUGSoutput$DIC

#==================================================================##
k <- 2
nMC <- 10000

res2_CF <- piv_MCMC(y = myLI2, k = k, nMC = nMC, 
                    software = "rjags")
rel2_CF <- piv_rel(res2_CF)

png(filename = "Parker_univ2_CF_trace_jags.png", width = 10, height = 4, units = "in",res=300)
piv_plot(y=myLI2, res2_CF, rel2_CF, par = c("mean"), type="chains")
dev.off()

png(filename = "Parker_univ2_CF_hist_jags.png", width = 10, height = 4, units = "in",res=300)
piv_plot(y=myLI2, res2_CF, rel2_CF, par = c("mean"), type="hist")
dev.off()


#==================================================================##
# 3 component Gaussian mixture model via Stan
#==================================================================##

#######
# Package the data for JAGS

data.package.CF3 <- list(N=length(myLI), y=myLI, k=3,
                         b0=0, B0inv=0.1,
                         nu0Half=1,
                         g0Half=1,g0G0Half=1, e=c(1,1,1))



params.to.monitor <- c("eta","mu","tau")

jags.fit.CF3 <- jags(data=data.package.CF3,inits=init.generator2,parameters.to.save=params.to.monitor,n.iter=10000,model.file="JAGS_mod1.txt",n.chains = 4,n.burnin = 2000,n.thin=5 )

jagsfitCF3.mcmc <- as.mcmc(jags.fit.CF3) 

summary(jagsfitCF3.mcmc)

plot(jagsfitCF3.mcmc[,"eta[1]"])
plot(jagsfitCF3.mcmc[,"eta[2]"])
plot(jagsfitCF3.mcmc[,"mu[1]"])
plot(jagsfitCF3.mcmc[,"mu[2]"])
plot(jagsfitCF3.mcmc[,"tau[1]"])
plot(jagsfitCF3.mcmc[,"tau[2]"])

densityplot(jagsfitCF3.mcmc)

DIC_CF3 <- jags.fit.CF3$BUGSoutput$DIC

#==================================================================##
k <- 3
nMC <- 10000

res3_CF <- piv_MCMC(y = myLI2, k = k, nMC = nMC, 
                    software = "rjags")
rel3_CF <- piv_rel(res3_CF)

png(filename = "Parker_univ3_CF_trace_jags.png", width = 10, height = 4, units = "in",res=300)
piv_plot(y=myLI2, res3_CF, rel3_CF, par = c("mean"), type="chains")
dev.off()

png(filename = "Parker_univ3_CF_hist_jags.png", width = 10, height = 4, units = "in",res=300)
piv_plot(y=myLI2, res3_CF, rel3_CF, par = c("mean"), type="hist")
dev.off()


#==================================================================##
# 4 component Gaussian mixture model via JAGS
#==================================================================##

#######
# Package the data for JAGS

data.package.CF4 <- list(N=length(myLI), y=myLI, k=4,
                         b0=0, B0inv=0.1,
                         nu0Half=1,
                         g0Half=1,g0G0Half=1, e=c(1,1,1,1))


params.to.monitor <- c("eta","mu","tau")

jags.fit.CF4 <- jags(data=data.package.CF4,inits=init.generator3,parameters.to.save=params.to.monitor,n.iter=10000,model.file="JAGS_mod1.txt",n.chains = 4,n.burnin = 2000,n.thin=5 )

jagsfitCF4.mcmc <- as.mcmc(jags.fit.CF4) 

summary(jagsfitCF4.mcmc)

plot(jagsfitCF4.mcmc[,"eta[1]"])
plot(jagsfitCF4.mcmc[,"eta[2]"])
plot(jagsfitCF4.mcmc[,"mu[1]"])
plot(jagsfitCF4.mcmc[,"mu[2]"])
plot(jagsfitCF4.mcmc[,"tau[1]"])
plot(jagsfitCF4.mcmc[,"tau[2]"])

densityplot(jagsfitCF4.mcmc)

DIC_CF4 <- jags.fit.CF4$BUGSoutput$DIC


#==================================================================##
k <- 4
nMC <- 10000

res4_CF <- piv_MCMC(y = myLI2, k = k, nMC = nMC, 
                    software = "rjags")
rel4_CF <- piv_rel(res4_CF)

png(filename = "Parker_univ4_CF_trace_jags.png", width = 10, height = 4, units = "in",res=300)
piv_plot(y=myLI2, res4_CF, rel4_CF, par = c("mean"), type="chains")
dev.off()

png(filename = "Parker_univ4_CF_hist_jags.png", width = 10, height = 4, units = "in",res=300)
piv_plot(y=myLI2, res4_CF, rel4_CF, par = c("mean"), type="hist")
dev.off()

save.image("/Volumes/PSYHOME/PSYRES/pthompson/DVMB/BGMM_laterality/mar03.Rdata")
savehistory("/Volumes/PSYHOME/PSYRES/pthompson/DVMB/BGMM_laterality/hist_mar03.Rhistory")

#==================================================================##
#COMPARE DICS

#DL

DIC_DL2
DIC_DL3
DIC_DL4


#CF

DIC_CF2
DIC_CF3
DIC_CF4