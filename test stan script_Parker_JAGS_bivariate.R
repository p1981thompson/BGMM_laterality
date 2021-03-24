library(pivmet)
library(bayesmix)
library(bayesplot)
library(tidyverse)
library(R2jags) 
library(runjags)
library(mixtools)
library(coda)
library(lattice)
library(MCMCvis)
library(MixtureInf)


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
myLI3<-LIs$Day1_finger_LI[complete.cases(LIs$Day1_finger_LI)]

#==================================================================##

LI_both_CF<-LIs[,c("Day1_CF_acc_LI","Day2_CF_acc_LI")] %>% remove_missing()

LI_both_DL<-LIs[,c("Day1_DL_acc_LI","Day2_DL_acc_LI")] %>% remove_missing()

LI_both_FT<-LIs[,c("Day1_finger_LI","Day2_finger_LI")] %>% remove_missing()


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

LI_fit_FT<-normalmixEM(myLI3)
plot(LI_fit_FT, density=TRUE, cex.axis=1.4, cex.lab=1.4, cex.main=1.8,main2="Finger tap", xlab2="LI")



LI.DL.boot <- boot.comp(myLI,max.comp=10,mix.type="normalmix", maxit=400,epsilon=1e-2)

LI.CF.boot <- boot.comp(myLI2,max.comp=10,mix.type="normalmix", maxit=400,epsilon=1e-2)

LI.FT.boot <- boot.comp(myLI3,max.comp=10,mix.type="normalmix", maxit=400,epsilon=1e-2)

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

#Finger Tapping#

LI.FT.k2 <- normalmixEM(myLI3,k=2,maxit=100,epsilon=0.01)
LI.FT.k3 <- normalmixEM(myLI3,k=3,maxit=100,epsilon=0.01)
LI.FT.k4 <- normalmixEM(myLI3,k=4,maxit=100,epsilon=0.01)

plot(hist(myLI3,breaks=101),col="grey",border="grey",freq=FALSE,
     xlab="LI",main="Finger Tapping (k=2)")
lines(density(myLI3),lty=2)
sapply(1:2,plot.normal.components,mixture=LI.FT.k2)


plot(hist(myLI3,breaks=101),col="grey",border="grey",freq=FALSE,
     xlab="LI",main="Finger Tapping (k=3)")
lines(density(myLI3),lty=2)
sapply(1:3,plot.normal.components,mixture=LI.FT.k3)

plot(hist(myLI3,breaks=101),col="grey",border="grey",freq=FALSE,
     xlab="LI",main="Finger Tapping (k=4)")
lines(density(myLI3),lty=2)
sapply(1:4,plot.normal.components,mixture=LI.FT.k4)


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
CVfun(myLI3)

#Pretty inconclusive!

MixtureInf::emtest.norm(myLI)
MixtureInf::emtest.norm(myLI2)
MixtureInf::emtest.norm(myLI3)

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
  
  mu[1] ~ dnorm(b01,B0inv);
  mu[2] ~ dnorm(b02,B0inv) T(mu[1],);
  
  for (j in 1:k) {
    
    tau[j] ~ dgamma(nu0Half,nu0S0Half);
    sigma[j] <- 1/tau[j]
  }
  S0 ~ dgamma(g0Half,g0G0Half);
  nu0S0Half <- nu0Half * S0;
  
  eta[1:k] ~ ddirch(e[1:k]);
}
",file="JAGS_mod1.txt")

cat("
model	{
  for (i in 1:N) {
    y[i] ~ dnorm(mu[S[i]],tau[S[i]]);
    S[i] ~ dcat(eta[]);
  }
  
  mu[1] ~ dnorm(b01,B0inv);
  mu[2] ~ dnorm(b02,B0inv) T(mu[1],);
  mu[3] ~ dnorm(b03,B0inv) T(mu[2],);
  
  for (j in 1:k) {
    
    tau[j] ~ dgamma(nu0Half,nu0S0Half);
    sigma[j] <- 1/tau[j]
  }
  S0 ~ dgamma(g0Half,g0G0Half);
  nu0S0Half <- nu0Half * S0;
  
  eta[1:k] ~ ddirch(e[1:k]);
}
",file="JAGS_mod2.txt")

cat("
model	{
  for (i in 1:N) {
    y[i] ~ dnorm(mu[S[i]],tau[S[i]]);
    S[i] ~ dcat(eta[]);
  }
  
  mu[1] ~ dnorm(b01,B0inv);
  mu[2] ~ dnorm(b02,B0inv) T(mu[1],);
  mu[3] ~ dnorm(b03,B0inv) T(mu[2],);
  mu[4] ~ dnorm(b04,B0inv) T(mu[3],);
  
  for (j in 1:k) {
    
    tau[j] ~ dgamma(nu0Half,nu0S0Half);
    sigma[j] <- 1/tau[j]
  }
  S0 ~ dgamma(g0Half,g0G0Half);
  nu0S0Half <- nu0Half * S0;
  
  eta[1:k] ~ ddirch(e[1:k]);
}
",file="JAGS_mod3.txt")


##########
# Make a function for generating initial guesses

init.generator1 <- function(){ list(
  mu = c(runif(1, -100,0),runif(1, 0,100)),
  tau = 1/runif(2, 1,10),
  eta = runif(2, 0,1))
}

init.generator2 <- function(){ list(
  mu = c(runif(1, -100,-33),runif(1, -32.9999,33),runif(1, 33.0000001,100)),
  tau = 1/runif(3, 1,10),
  eta = runif(3, 0,1))
}

init.generator3 <- function(){ list(
  mu = c(runif(1, -100,-50),runif(1, -49.9999,0),runif(1, 0.00001,50),runif(1, 50.00001,100)),
  tau = 1/runif(4, 1,10),
  eta = runif(4, 0,1))
}

#==================================================================##
# Dichotic listening
#==================================================================##
# 2 component Gaussian mixture model via jags
#==================================================================##
#
#######
# Package the data for JAGS

data.package.DL2 <- list(N=length(myLI), y=myLI, k=2,
                         b01=-20,b02=20, B0inv=0.05,
                         nu0Half=1,
                         g0Half=1,g0G0Half=20, e=c(1,1))


params.to.monitor <- c("eta","mu","tau","S")

jags.fit.DL2 <- jags(data=data.package.DL2,inits=init.generator1,parameters.to.save=params.to.monitor,n.iter=20000,model.file="JAGS_mod1.txt",n.chains = 4,n.burnin = 4000,n.thin=5 )

jagsfitDL2.mcmc <- as.mcmc(jags.fit.DL2) 

summary(jagsfitDL2.mcmc)

plot(jagsfitDL2.mcmc[,"eta[1]"])
plot(jagsfitDL2.mcmc[,"eta[2]"])
plot(jagsfitDL2.mcmc[,"mu[1]"])
plot(jagsfitDL2.mcmc[,"mu[2]"])
plot(jagsfitDL2.mcmc[,"tau[1]"])
plot(jagsfitDL2.mcmc[,"tau[2]"])

#densityplot(jagsfitDL2.mcmc,parameters=c("eta","mu","tau"))

DIC_DL2 <- jags.fit.DL2$BUGSoutput$DIC



#==================================================================##
k <- 2
nMC <- 20000

res2_DL <- piv_MCMC(y = myLI, k = k, nMC = nMC, software = "rjags", priors=list(b0=0, B0inv=0.05, nu0Half=1, g0Half=1, g0G0Half=20, e=c(1,1)))

rel2_DL <- piv_rel(res2_DL)

png(filename = "Parker_univ2_DL_trace_jags.png", width = 10, height = 4, units = "in",res=300)
piv_plot(y=myLI, res2_DL, rel2_DL, par = c("mean"), type="chains")
dev.off()

png(filename = "Parker_univ2_DL_hist_jags.png", width = 10, height = 4, units = "in",res=300)
piv_plot(y=myLI, res2_DL, rel2_DL, par = c("mean"), type="hist")
dev.off()

#==================================================================##
# 3 component Gaussian mixture model via jags
#==================================================================##

#######
# Package the data for JAGS

data.package.DL3 <- list(N=length(myLI), y=myLI, k=3,
                         b01=-50, b02=0, b03=50, B0inv=0.05,
                         nu0Half=1,
                         g0Half=1,g0G0Half=20, e=c(1,1,1))

params.to.monitor <- c("eta","mu","tau")

jags.fit.DL3 <- jags(data=data.package.DL3,inits=init.generator2,parameters.to.save=params.to.monitor,n.iter=100000,model.file="JAGS_mod2.txt",n.chains = 4,n.burnin = 20000,n.thin=5 )

jagsfitDL3.mcmc <- as.mcmc(jags.fit.DL3) 

summary(jagsfitDL3.mcmc)

plot(jagsfitDL3.mcmc[,"eta[1]"])
plot(jagsfitDL3.mcmc[,"eta[2]"])
plot(jagsfitDL3.mcmc[,"eta[3]"])
plot(jagsfitDL3.mcmc[,"mu[1]"])
plot(jagsfitDL3.mcmc[,"mu[2]"])
plot(jagsfitDL3.mcmc[,"mu[3]"])
plot(jagsfitDL3.mcmc[,"tau[1]"])
plot(jagsfitDL3.mcmc[,"tau[2]"])
plot(jagsfitDL3.mcmc[,"tau[3]"])

#densityplot(jagsfitDL3.mcmc,parameters=c("eta","mu","tau"))

DIC_DL3 <- jags.fit.DL3$BUGSoutput$DIC


#==================================================================##
k <- 3
nMC <- 20000

res3_DL <- piv_MCMC(y = myLI, k = k, nMC = nMC, 
                    software = "rjags",priors=list(b0=0, B0inv=0.05,
                                                   nu0Half=1,
                                                   g0Half=1,g0G0Half=20, e=c(1,1,1)))
rel3_DL <- piv_rel(res3_DL)

png(filename = "Parker_univ3_DL_trace_jags.png", width = 10, height = 4, units = "in",res=300)
piv_plot(y=myLI, res3_DL, rel3_DL, par = c("mean"), type="chains")
dev.off()

png(filename = "Parker_univ3_DL_hist_jags.png", width = 10, height = 4, units = "in",res=300)
piv_plot(y=myLI, res3_DL, rel3_DL, par = c("mean"), type="hist")
dev.off()

#==================================================================##
# 4 component Gaussian mixture model via jags
#==================================================================##

#######
# Package the data for JAGS

data.package.DL4 <- list(N=length(myLI), y=myLI, k=4,
                         b01=-20,b02=-10,b03=10,b04=20, B0inv=0.05,
                         nu0Half=1,
                         g0Half=1,g0G0Half=20, e=c(1,1,1,1))

params.to.monitor <- c("eta","mu","tau")

jags.fit.DL4 <- jags(data=data.package.DL4,inits=init.generator3,parameters.to.save=params.to.monitor,n.iter=20000,model.file="JAGS_mod3.txt",n.chains = 4,n.burnin = 4000,n.thin=5 )

jagsfitDL4.mcmc <- as.mcmc(jags.fit.DL4) 

summary(jagsfitDL4.mcmc)

plot(jagsfitDL4.mcmc[,"eta[1]"])
plot(jagsfitDL4.mcmc[,"eta[2]"])
plot(jagsfitDL4.mcmc[,"eta[3]"])
plot(jagsfitDL4.mcmc[,"eta[4]"])
plot(jagsfitDL4.mcmc[,"mu[1]"])
plot(jagsfitDL4.mcmc[,"mu[2]"])
plot(jagsfitDL4.mcmc[,"mu[3]"])
plot(jagsfitDL4.mcmc[,"mu[4]"])
plot(jagsfitDL4.mcmc[,"tau[1]"])
plot(jagsfitDL4.mcmc[,"tau[2]"])
plot(jagsfitDL4.mcmc[,"tau[3]"])
plot(jagsfitDL4.mcmc[,"tau[4]"])

#densityplot(jagsfitDL4.mcmc,parameters=c("eta","mu","tau"))

DIC_DL4 <- jags.fit.DL4$BUGSoutput$DIC

#==================================================================##
k <- 4
nMC <- 20000

res4_DL <- piv_MCMC(y = myLI, k = k, nMC = nMC, 
                    software = "rjags",priors=list(b0=0, B0inv=0.05,
                                                   nu0Half=1,
                                                   g0Half=1,g0G0Half=20, e=c(1,1,1,1)))
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
# 2 component Gaussian mixture model via jags
#==================================================================##

#######
# Package the data for JAGS

data.package.CF2 <- list(N=length(myLI), y=myLI, k=2,
                         b01=-50,b02=50, B0inv=0.05,
                         nu0Half=1,
                         g0Half=1,g0G0Half=20, e=c(1,1))

params.to.monitor <- c("eta","mu","tau")

jags.fit.CF2 <- jags(data=data.package.CF2,inits=init.generator1,parameters.to.save=params.to.monitor,n.iter=20000,model.file="JAGS_mod1.txt",n.chains = 4,n.burnin = 4000,n.thin=5 )

jagsfitCF2.mcmc <- as.mcmc(jags.fit.CF2) 

summary(jagsfitCF2.mcmc)

plot(jagsfitCF2.mcmc[,"eta[1]"])
plot(jagsfitCF2.mcmc[,"eta[2]"])
plot(jagsfitCF2.mcmc[,"mu[1]"])
plot(jagsfitCF2.mcmc[,"mu[2]"])
plot(jagsfitCF2.mcmc[,"tau[1]"])
plot(jagsfitCF2.mcmc[,"tau[2]"])

#densityplot(jagsfitCF2.mcmc,parameters=c("eta","mu","tau"))

DIC_CF2 <- jags.fit.CF2$BUGSoutput$DIC

#==================================================================##
k <- 2
nMC <- 20000

res2_CF <- piv_MCMC(y = myLI2, k = k, nMC = nMC, 
                    software = "rjags",priors=list(b0=0, B0inv=0.05,
                                                   nu0Half=1,
                                                   g0Half=1,g0G0Half=20, e=c(1,1)))
rel2_CF <- piv_rel(res2_CF)

png(filename = "Parker_univ2_CF_trace_jags.png", width = 10, height = 4, units = "in",res=300)
piv_plot(y=myLI2, res2_CF, rel2_CF, par = c("mean"), type="chains")
dev.off()

png(filename = "Parker_univ2_CF_hist_jags.png", width = 10, height = 4, units = "in",res=300)
piv_plot(y=myLI2, res2_CF, rel2_CF, par = c("mean"), type="hist")
dev.off()


#==================================================================##
# 3 component Gaussian mixture model via jags
#==================================================================##

#######
# Package the data for JAGS

data.package.CF3 <- list(N=length(myLI), y=myLI, k=3,
                         b01=-50,b02=0,b03=50,B0inv=0.05,
                         nu0Half=1,
                         g0Half=1,g0G0Half=20, e=c(1,1,1))



params.to.monitor <- c("eta","mu","tau")

jags.fit.CF3 <- jags(data=data.package.CF3,inits=init.generator2,parameters.to.save=params.to.monitor,n.iter=20000,model.file="JAGS_mod2.txt",n.chains = 4,n.burnin = 4000,n.thin=5 )

jagsfitCF3.mcmc <- as.mcmc(jags.fit.CF3) 

summary(jagsfitCF3.mcmc)

plot(jagsfitCF3.mcmc[,"eta[1]"])
plot(jagsfitCF3.mcmc[,"eta[2]"])
plot(jagsfitCF3.mcmc[,"mu[1]"])
plot(jagsfitCF3.mcmc[,"mu[2]"])
plot(jagsfitCF3.mcmc[,"tau[1]"])
plot(jagsfitCF3.mcmc[,"tau[2]"])

#densityplot(jagsfitCF3.mcmc,parameters=c("eta","mu","tau"))

DIC_CF3 <- jags.fit.CF3$BUGSoutput$DIC

#==================================================================##
k <- 3
nMC <- 20000

res3_CF <- piv_MCMC(y = myLI2, k = k, nMC = nMC, software = "rjags",priors=list(b0=0, B0inv=0.05,nu0Half=1,g0Half=1,g0G0Half=20, e=c(1,1,1)))
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
                         b01=-50,b02=-20,b03=20,b04=50,B0inv=0.05,
                         nu0Half=1,
                         g0Half=1,g0G0Half=20, e=c(1,1,1,1))


params.to.monitor <- c("eta","mu","tau")

jags.fit.CF4 <- jags(data=data.package.CF4,inits=init.generator3,parameters.to.save=params.to.monitor,n.iter=20000,model.file="JAGS_mod3.txt",n.chains = 4,n.burnin = 4000,n.thin=5 )

jagsfitCF4.mcmc <- as.mcmc(jags.fit.CF4) 

summary(jagsfitCF4.mcmc)

plot(jagsfitCF4.mcmc[,"eta[1]"])
plot(jagsfitCF4.mcmc[,"eta[2]"])
plot(jagsfitCF4.mcmc[,"mu[1]"])
plot(jagsfitCF4.mcmc[,"mu[2]"])
plot(jagsfitCF4.mcmc[,"tau[1]"])
plot(jagsfitCF4.mcmc[,"tau[2]"])

#densityplot(jagsfitCF4.mcmc,parameters=c("eta","mu","tau"))

DIC_CF4 <- jags.fit.CF4$BUGSoutput$DIC


#==================================================================##
k <- 4
nMC <- 20000

res4_CF <- piv_MCMC(y = myLI2, k = k, nMC = nMC, software = "rjags",priors=list(b0=0, B0inv=0.05,nu0Half=1,g0Half=1,g0G0Half=20, e=c(1,1,1,1)))
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


#==================================================================###==================================================================##
#==================================================================###==================================================================##
#==================================================================###==================================================================##

# Bivariate mixture model for both time points combined.

# 08-03-2021
# 19-03-2021: Update.
#==================================================================###==================================================================###
#==================================================================###==================================================================###
#==================================================================###==================================================================###

# Model script for JAGS (2 class model)  

cat(
  "model{
# Likelihood:
for (i in 1:N){
yprev[i,1:D]<-y[i,1:D]
y[i,1:D] ~ dmnorm(mu[clust[i],],tau)
clust[i] ~ dcat(eta[1:k] )
}

# Prior:
mu[1,1:D] ~ dmnorm(mu_0A[],S2[,])
mu[2,1:D] ~ dmnorm(mu_0B[],S2[,])
tau[1:D,1:D] ~ dwish(S3[,],D+1)
Sigma[1:D,1:D] <- inverse(tau[,])
eta[1:k] ~ ddirch(alpha)  
}",
file="JAGS_bivar.txt"  
)
#==================================================================##

# Model script for JAGS (2 class model)  

cat(
  "model{
# Likelihood:
for (i in 1:N){
yprev[i,1:D]<-y[i,1:D]
y[i,1:D] ~ dmnorm(mu[clust[i],],tau)
clust[i] ~ dcat(eta[1:k] )
}

# Prior:
mu[1,1:D] ~ dmnorm(mu_0A[],S2[,])
mu[2,1:D] ~ dmnorm(mu_0B[],S2[,])
mu[3,1:D] ~ dmnorm(mu_0C[],S2[,])
tau[1:D,1:D] ~ dwish(S3[,],D+1)
Sigma[1:D,1:D] <- inverse(tau[,])
eta[1:k] ~ ddirch(alpha)  
}",
file="JAGS_bivar3.txt"  
)
#==================================================================##

# Model script for JAGS (2 class model)  

cat(
  "model{
# Likelihood:
for (i in 1:N){
yprev[i,1:D]<-y[i,1:D]
y[i,1:D] ~ dmnorm(mu[clust[i],],tau)
clust[i] ~ dcat(eta[1:k] )
}

# Prior:
mu[1,1:D] ~ dmnorm(mu_0A[],S2[,])
mu[2,1:D] ~ dmnorm(mu_0B[],S2[,])
mu[3,1:D] ~ dmnorm(mu_0C[],S2[,])
mu[4,1:D] ~ dmnorm(mu_0D[],S2[,])
tau[1:D,1:D] ~ dwish(S3[,],D+1)
Sigma[1:D,1:D] <- inverse(tau[,])
eta[1:k] ~ ddirch(alpha)  
}",
file="JAGS_bivar4.txt"  
)

#==================================================================##

#Generates initial values for the JAGS MCMC sampler. (this function is specific to the two class model)

init.generator <- function(){ list(
  mu = matrix(c(runif(1, -100,0),runif(1, 0,100),runif(1, -100,0),runif(1, 0,100)),2,2,byrow=TRUE),
  tau = matrix(c(1/runif(1, 1,10),0,0,1/runif(1, 1,10)),2,2,byrow=TRUE),
  eta = runif(2, 0,1))
}

#==================================================================##
# Script components for running the JAGS model.
k=2 # No. of classes
D=2 # No. of time points
#==================================================================##

# Set up the data objects for the JAGS models. This includes the data and fixed parameters (i.e. N, k, and D). It also has information on hyper-parameters within the priors. 
data.package.CF2x2 <- list(N=dim(LI_both_CF)[1],D=dim(LI_both_CF)[2], y=LI_both_CF, k=2,
                           S2= diag(D)*1000, S3= diag(D)*1000, mu_0A=rep(-50, D), mu_0B=rep(50, D), alpha = rep(1,k))

# selects which parameters to follow: eta - mixing weights, mu - class means, tau - 1/variance per class, clust - vector of class assignments of individual 
params.to.monitor <- c("eta","mu","tau","clust")

# Fits the two class, bivariate mixture model (bivariate as both time points per measure).
jags.fit.CF2x2 <- jags(data=data.package.CF2x2,inits=init.generator,parameters.to.save=params.to.monitor,n.iter=20000,model.file="JAGS_bivar.txt",n.chains = 4,n.burnin = 4000,n.thin=5 )

# Extracts MCMC chains
jagsfitCF2x2.mcmc <- as.mcmc(jags.fit.CF2x2) 

summary(jagsfitCF2x2.mcmc)

# MCMC trace plots for checking stability of posterior sampling.
plot(jagsfitCF2x2.mcmc[,"eta[1]"])
plot(jagsfitCF2x2.mcmc[,"eta[2]"])
plot(jagsfitCF2x2.mcmc[,"mu[1,1]"])
plot(jagsfitCF2x2.mcmc[,"mu[1,2]"])
plot(jagsfitCF2x2.mcmc[,"mu[2,1]"])
plot(jagsfitCF2x2.mcmc[,"mu[2,2]"])
plot(jagsfitCF2x2.mcmc[,"tau[1,1]"])
plot(jagsfitCF2x2.mcmc[,"tau[2,1]"])
plot(jagsfitCF2x2.mcmc[,"tau[1,2]"])
plot(jagsfitCF2x2.mcmc[,"tau[2,2]"])

#densityplot(jagsfitDL2.mcmc,parameters=c("eta","mu","tau"))

# Deviance information criterion for model comparisons.
DIC_CF2x2 <- jags.fit.CF2x2$BUGSoutput$DIC

#==================================================================##
# Dichotic Listening (2 class, 2 time points)
#==================================================================##

#set up the data components including dataframe, and fixed parameters (i.e. N, D, k and hyperparameters.)
data.package.DL2x2 <- list(N=dim(LI_both_DL)[1],D=dim(LI_both_DL)[2], y=LI_both_DL, k=2,
                           S2= diag(D)*1000, S3= diag(D)*1000, mu_0A=rep(-50, D), mu_0B=rep(50, D), alpha = rep(1,k))

# Parameters to model
params.to.monitor <- c("eta","mu","tau","clust")

# Runs the JAGS model for dichotic listening, bivariate mixture model (2 class, 2 time points)
jags.fit.DL2x2 <- jags(data=data.package.DL2x2,inits=init.generator,parameters.to.save=params.to.monitor,n.iter=20000,model.file="JAGS_bivar.txt",n.chains = 4,n.burnin = 4000,n.thin=5 )

# extract MCMC samples for posterior evaluation
jagsfitDL2x2.mcmc <- as.mcmc(jags.fit.DL2x2) 

summary(jagsfitDL2x2.mcmc)

# plot MCMC chains for parameters monitored to check stability
plot(jagsfitDL2x2.mcmc[,"eta[1]"])
plot(jagsfitDL2x2.mcmc[,"eta[2]"])
plot(jagsfitDL2x2.mcmc[,"mu[1,1]"])
plot(jagsfitDL2x2.mcmc[,"mu[1,2]"])
plot(jagsfitDL2x2.mcmc[,"mu[2,1]"])
plot(jagsfitDL2x2.mcmc[,"mu[2,2]"])
plot(jagsfitDL2x2.mcmc[,"tau[1,1]"])
plot(jagsfitDL2x2.mcmc[,"tau[2,1]"])
plot(jagsfitDL2x2.mcmc[,"tau[1,2]"])
plot(jagsfitDL2x2.mcmc[,"tau[2,2]"])

#densityplot(jagsfitDL2.mcmc,parameters=c("eta","mu","tau"))

# Deviance information criterion for model comparison
DIC_DL2x2 <- jags.fit.DL2x2$BUGSoutput$DIC

#==================================================================##
# Finger Tapping (2 class, 2 time points)
#==================================================================##

# set up data components list (includes data and fixed parameters)
data.package.FT2x2 <- list(N=dim(LI_both_FT)[1],D=dim(LI_both_FT)[2], y=LI_both_FT, k=2,
                           S2= diag(D)*1000, S3= diag(D)*1000, mu_0A=rep(-50, D), mu_0B=rep(50, D), alpha = rep(1,k))

# Parameters to monitor in the output.
params.to.monitor <- c("eta","mu","tau","clust")

# Fit the JAGS model to the data (bivariate mixture model for finger tapping [2 class, 2 time points])
jags.fit.FT2x2 <- jags(data=data.package.FT2x2,inits=init.generator,parameters.to.save=params.to.monitor,n.iter=20000,model.file="JAGS_bivar.txt",n.chains = 4,n.burnin = 4000,n.thin=5 )

# extract MCMC samples  
jagsfitFT2x2.mcmc <- as.mcmc(jags.fit.DL2x2) 

summary(jagsfitFT2x2.mcmc)

# MCMC trace plots to check stability
plot(jagsfitFT2x2.mcmc[,"eta[1]"])
plot(jagsfitFT2x2.mcmc[,"eta[2]"])
plot(jagsfitFT2x2.mcmc[,"mu[1,1]"])
plot(jagsfitFT2x2.mcmc[,"mu[1,2]"])
plot(jagsfitFT2x2.mcmc[,"mu[2,1]"])
plot(jagsfitFT2x2.mcmc[,"mu[2,2]"])
plot(jagsfitFT2x2.mcmc[,"tau[1,1]"])
plot(jagsfitFT2x2.mcmc[,"tau[2,1]"])
plot(jagsfitFT2x2.mcmc[,"tau[1,2]"])
plot(jagsfitFT2x2.mcmc[,"tau[2,2]"])

#densityplot(jagsfitDL2.mcmc,parameters=c("eta","mu","tau"))

# Extract deviance information criterion.
DIC_FT2x2 <- jags.fit.FT2x2$BUGSoutput$DIC


#==================================================================##
# Repeat analysis using pivmet package. This ensures label switching 
# is fixed. The analysis is virtually identical model but the previous
# model above uses starting values and more specific priors.
 
#==================================================================##
# Chimeric faces
#==================================================================##

# run model with label switching adjustment
res_CF_both2<-piv_MCMC(y = as.matrix(LI_both_CF), k = 2, nMC = 10000, software = "rjags")
rel_CF_both2 <- piv_rel(mcmc=res_CF_both2) #extract relative estimates

piv_plot(y = as.matrix(LI_both_CF), mcmc = res_CF_both2, rel_est = rel_CF_both2, par = "mean", type = "chains")
piv_plot(y = as.matrix(LI_both_CF), mcmc = res_CF_both2, rel_est = rel_CF_both2, par = "weight", type = "chains")
#> Description: traceplot of the raw MCMC chains and the relabelled chains for the means parameters (coordinate 1 and 2). Each colored chain corresponds to one of the k distinct parameters of the mixture model. Overlapping chains may reveal that the MCMC sample is not able to distinguish between the components.
piv_plot(y = as.matrix(LI_both_CF), mcmc = res_CF_both2, rel_est = rel_CF_both2, type = "hist")
#> Description: 3d histogram of the data along with the posterior estimates of the relabelled means (triangle points)

#==================================================================##
#Dichotic Listening
#==================================================================##
# run model with label switching adjustment
res_DL_both2<-piv_MCMC(y = as.matrix(LI_both_DL), k = 2, nMC = 10000, software = "rjags")
rel_DL_both2 <- piv_rel(mcmc=res_DL_both2)#extract relative estimates

piv_plot(y = as.matrix(LI_both_DL), mcmc = res_DL_both2, rel_est = rel_DL_both2, par = "mean", type = "chains")
piv_plot(y = as.matrix(LI_both_DL), mcmc = res_DL_both2, rel_est = rel_DL_both2, par = "weight", type = "chains")
#> Description: traceplot of the raw MCMC chains and the relabelled chains for the means parameters (coordinate 1 and 2). Each colored chain corresponds to one of the k distinct parameters of the mixture model. Overlapping chains may reveal that the MCMC sample is not able to distinguish between the components.
piv_plot(y = as.matrix(LI_both_DL), mcmc = res_DL_both2, rel_est = rel_DL_both2, type = "hist")
#> Description: 3d histogram of the data along with the posterior estimates of the relabelled means (triangle points)

#==================================================================##
#Finger tapping
#==================================================================##
# run model with label switching adjustment
res_FT_both2<-piv_MCMC(y = as.matrix(LI_both_FT), k = 2, nMC = 10000, software = "rjags")
rel_FT_both2 <- piv_rel(mcmc=res_FT_both2)#extract relative estimates

piv_plot(y = as.matrix(LI_both_FT), mcmc = res_FT_both2, rel_est = rel_FT_both2, par = "mean", type = "chains")
piv_plot(y = as.matrix(LI_both_FT), mcmc = res_FT_both2, rel_est = rel_FT_both2, par = "weight", type = "chains")
#> Description: traceplot of the raw MCMC chains and the relabelled chains for the means parameters (coordinate 1 and 2). Each colored chain corresponds to one of the k distinct parameters of the mixture model. Overlapping chains may reveal that the MCMC sample is not able to distinguish between the components.
piv_plot(y = as.matrix(LI_both_FT), mcmc = res_FT_both2, rel_est = rel_FT_both2, type = "hist")
#> Description: 3d histogram of the data along with the posterior estimates of the relabelled means (triangle points)
#==================================================================##
#==================================================================##
#==================================================================##

# THREE CLASS MODELS!

#==================================================================##
#######
# Package the data for JAGS
k=3 # 3-classes
D=2 # two time points

# New function to generate initial values based on a 3-class model
init.generator <- function(){ list(
  mu = matrix(c(runif(1, -100,-33),runif(1, -33.000001,33),runif(1, 33.00001,100),runif(1, -100,-33),runif(1, -33.000001,33),runif(1, 33.00001,100)),3,2,byrow=F),
  tau = matrix(c(1/runif(1, 1,10),0,0,1/runif(1, 1,10)),2,2,byrow=TRUE),
  eta = runif(3, 0,1))
}
#==================================================================##

# creates the list of data components ( data, fixed parameters and hyperparameters in priors.)
data.package.CF2x3 <- list(N=dim(LI_both_CF)[1],D=dim(LI_both_CF)[2], y=LI_both_CF, k=3,
                           S2= diag(D)*1000, S3= diag(D)*1000, mu_0A=rep(-50, D), mu_0B=rep(0, D),mu_0C=rep(50, D), alpha = rep(1,k))

# parameters to follow.
params.to.monitor <- c("eta","mu","tau","clust")

# fits the JAGS bivariate mixture model for 3-class chimeric faces. 
jags.fit.CF2x3 <- jags(data=data.package.CF2x3,inits=init.generator,parameters.to.save=params.to.monitor,n.iter=20000,model.file="JAGS_bivar3.txt",n.chains = 4,n.burnin = 4000,n.thin=5 )

# Extract MCMC samples.
jagsfitCF2x3.mcmc <- as.mcmc(jags.fit.CF2x3) 

summary(jagsfitCF2x3.mcmc)

# MCMC trace plots for the 3 class CF model parameters 
plot(jagsfitCF2x3.mcmc[,"eta[1]"])
plot(jagsfitCF2x3.mcmc[,"eta[2]"])
plot(jagsfitCF2x3.mcmc[,"eta[3]"])
plot(jagsfitCF2x3.mcmc[,"mu[1,1]"])
plot(jagsfitCF2x3.mcmc[,"mu[1,2]"])
plot(jagsfitCF2x3.mcmc[,"mu[2,1]"])
plot(jagsfitCF2x3.mcmc[,"mu[2,2]"])
plot(jagsfitCF2x3.mcmc[,"mu[3,1]"])
plot(jagsfitCF2x3.mcmc[,"mu[3,2]"])
plot(jagsfitCF2x3.mcmc[,"tau[1,1]"])
plot(jagsfitCF2x3.mcmc[,"tau[2,1]"])
plot(jagsfitCF2x3.mcmc[,"tau[1,2]"])
plot(jagsfitCF2x3.mcmc[,"tau[2,2]"])

#densityplot(jagsfitDL2.mcmc,parameters=c("eta","mu","tau"))

# Deviance Information Criterion for model comparison
DIC_CF2x3 <- jags.fit.CF2x3$BUGSoutput$DIC

#==================================================================##
# Dichotic Listening (3 class, 2 time points)
#==================================================================##
# creates the list of data components ( data, fixed parameters and hyperparameters in priors.)
data.package.DL2x3 <- list(N=dim(LI_both_DL)[1],D=dim(LI_both_DL)[2], y=LI_both_DL, k=3,
                           S2= diag(D)*1000, S3= diag(D)*1000, mu_0A=rep(-50, D), mu_0B=rep(0, D), mu_0C=rep(50, D), alpha = rep(1,k))

# parameters to follow.
params.to.monitor <- c("eta","mu","tau","clust")

# fits the JAGS bivariate mixture model for 3-class Dichotic listening
jags.fit.DL2x3 <- jags(data=data.package.DL2x3,inits=init.generator,parameters.to.save=params.to.monitor,n.iter=20000,model.file="JAGS_bivar3.txt",n.chains = 4,n.burnin = 4000,n.thin=5 )

# Extract MCMC samples.
jagsfitDL2x3.mcmc <- as.mcmc(jags.fit.DL2x3) 

summary(jagsfitDL2x3.mcmc)

# MCMC trace plots for the 3 class DL model parameters 
plot(jagsfitDL2x3.mcmc[,"eta[1]"])
plot(jagsfitDL2x3.mcmc[,"eta[2]"])
plot(jagsfitDL2x3.mcmc[,"eta[3]"])
plot(jagsfitDL2x3.mcmc[,"mu[1,1]"])
plot(jagsfitDL2x3.mcmc[,"mu[1,2]"])
plot(jagsfitDL2x3.mcmc[,"mu[2,1]"])
plot(jagsfitDL2x3.mcmc[,"mu[2,2]"])
plot(jagsfitDL2x3.mcmc[,"mu[3,1]"])
plot(jagsfitDL2x3.mcmc[,"mu[3,2]"])
plot(jagsfitDL2x3.mcmc[,"tau[1,1]"])
plot(jagsfitDL2x3.mcmc[,"tau[2,1]"])
plot(jagsfitDL2x3.mcmc[,"tau[1,2]"])
plot(jagsfitDL2x3.mcmc[,"tau[2,2]"])

#densityplot(jagsfitDL2.mcmc,parameters=c("eta","mu","tau"))

# Deviance Information Criterion for model comparison
DIC_DL2x3 <- jags.fit.DL2x3$BUGSoutput$DIC

#==================================================================##
# Finger Tapping (2 class, 2 time points)
#==================================================================##
# creates the list of data components ( data, fixed parameters and hyperparameters in priors.)
data.package.FT2x3 <- list(N=dim(LI_both_FT)[1],D=dim(LI_both_FT)[2], y=LI_both_FT, k=3,
                           S2= diag(D)*1000, S3= diag(D)*1000, mu_0A=rep(-50, D), mu_0B=rep(0, D), mu_0C=rep(50, D), alpha = rep(1,k))

# parameters to follow.
params.to.monitor <- c("eta","mu","tau","clust")

# fits the JAGS bivariate mixture model for 3-class Finger tapping
jags.fit.FT2x3 <- jags(data=data.package.FT2x3,inits=init.generator,parameters.to.save=params.to.monitor,n.iter=20000,model.file="JAGS_bivar3.txt",n.chains = 4,n.burnin = 4000,n.thin=5 )

# Extract MCMC samples.
jagsfitFT2x3.mcmc <- as.mcmc(jags.fit.DL3x2) 

summary(jagsfitFT2x3.mcmc)

# MCMC trace plots for the 3 class FT model parameters 
plot(jagsfitFT2x3.mcmc[,"eta[1]"])
plot(jagsfitFT2x3.mcmc[,"eta[2]"])
plot(jagsfitFT2x3.mcmc[,"eta[3]"])
plot(jagsfitFT2x3.mcmc[,"mu[1,1]"])
plot(jagsfitFT2x3.mcmc[,"mu[1,2]"])
plot(jagsfitFT2x3.mcmc[,"mu[2,1]"])
plot(jagsfitFT2x3.mcmc[,"mu[2,2]"])
plot(jagsfitFT2x3.mcmc[,"mu[3,1]"])
plot(jagsfitFT2x3.mcmc[,"mu[3,2]"])
plot(jagsfitFT2x3.mcmc[,"tau[1,1]"])
plot(jagsfitFT2x3.mcmc[,"tau[2,1]"])
plot(jagsfitFT2x3.mcmc[,"tau[1,2]"])
plot(jagsfitFT2x3.mcmc[,"tau[2,2]"])

#densityplot(jagsfitDL2.mcmc,parameters=c("eta","mu","tau"))

# Deviance Information Criterion for model comparison
DIC_FT2x3 <- jags.fit.FT2x3$BUGSoutput$DIC

#==================================================================##
# Repeat analysis using pivmet package. This ensures label switching 
# is fixed. The analysis is virtually identical model but the previous
# model above uses starting values and more specific priors.

#==================================================================##
# Chimeric faces
#==================================================================##

# run model with label switching adjustment
res_CF_both3<-piv_MCMC(y = as.matrix(LI_both_CF), k = 3, nMC = 10000, software = "rjags")
rel_CF_both3 <- piv_rel(mcmc=res_CF_both3) #extract relative estimates

piv_plot(y = as.matrix(LI_both_CF), mcmc = res_CF_both3, rel_est = rel_CF_both3, par = "mean", type = "chains")
piv_plot(y = as.matrix(LI_both_CF), mcmc = res_CF_both3, rel_est = rel_CF_both3, par = "weight", type = "chains")
#> Description: traceplot of the raw MCMC chains and the relabelled chains for the means parameters (coordinate 1 and 2). Each colored chain corresponds to one of the k distinct parameters of the mixture model. Overlapping chains may reveal that the MCMC sample is not able to distinguish between the components.
piv_plot(y = as.matrix(LI_both_CF), mcmc = res_CF_both3, rel_est = rel_CF_both3, type = "hist")
#> Description: 3d histogram of the data along with the posterior estimates of the relabelled means (triangle points)

#==================================================================##
#Dichotic Listening
#==================================================================##
# run model with label switching adjustment
res_DL_both3<-piv_MCMC(y = as.matrix(LI_both_DL), k = 3, nMC = 10000, software = "rjags")
rel_DL_both3 <- piv_rel(mcmc=res_DL_both3)#extract relative estimates

piv_plot(y = as.matrix(LI_both_DL), mcmc = res_DL_both3, rel_est = rel_DL_both3, par = "mean", type = "chains")
piv_plot(y = as.matrix(LI_both_DL), mcmc = res_DL_both3, rel_est = rel_DL_both3, par = "weight", type = "chains")
#> Description: traceplot of the raw MCMC chains and the relabelled chains for the means parameters (coordinate 1 and 2). Each colored chain corresponds to one of the k distinct parameters of the mixture model. Overlapping chains may reveal that the MCMC sample is not able to distinguish between the components.
piv_plot(y = as.matrix(LI_both_DL), mcmc = res_DL_both3, rel_est = rel_DL_both3, type = "hist")
#> Description: 3d histogram of the data along with the posterior estimates of the relabelled means (triangle points)

#==================================================================##
#Finger tapping
#==================================================================##
# run model with label switching adjustment
res_FT_both3<-piv_MCMC(y = as.matrix(LI_both_FT), k = 3, nMC = 10000, software = "rjags")
rel_FT_both3 <- piv_rel(mcmc=res_FT_both3)#extract relative estimates

piv_plot(y = as.matrix(LI_both_FT), mcmc = res_FT_both3, rel_est = rel_FT_both3, par = "mean", type = "chains")
piv_plot(y = as.matrix(LI_both_FT), mcmc = res_FT_both3, rel_est = rel_FT_both3, par = "weight", type = "chains")
#> Description: traceplot of the raw MCMC chains and the relabelled chains for the means parameters (coordinate 1 and 2). Each colored chain corresponds to one of the k distinct parameters of the mixture model. Overlapping chains may reveal that the MCMC sample is not able to distinguish between the components.
piv_plot(y = as.matrix(LI_both_FT), mcmc = res_FT_both3, rel_est = rel_FT_both3, type = "hist")
#> Description: 3d histogram of the data along with the posterior estimates of the relabelled means (triangle points)
#==================================================================##
#==================================================================##
#==================================================================##

# THREE CLASS MODELS!

#==================================================================##
# New function to generate initial values based on a 3-class model
# Package the data for JAGS (4 class [k=4], both time points [D=2])
k=4
D=2

# generate starting values for the sampler
init.generator <- function(){ list(
  mu = matrix(c(runif(1, -100,-50),runif(1, -50.000001,0),runif(1, 0.00001,50),runif(1, 50,100),runif(1, -100,-50),runif(1, -50.000001,0),runif(1, 0.00001,50),runif(1, 50,100)),4,2,byrow=F),
  tau = matrix(c(1/runif(1, 1,10),0,0,1/runif(1, 1,10)),2,2,byrow=TRUE),
  eta = runif(4, 0,1))
}
#==================================================================##
#set up data
data.package.CF2x4 <- list(N=dim(LI_both_CF)[1],D=dim(LI_both_CF)[2], y=LI_both_CF, k=4,
                           S2= diag(D)*1000, S3= diag(D)*1000, mu_0A=rep(-75, D), mu_0B=rep(-25, D),mu_0C=rep(25, D),mu_0D=rep(75, D), alpha = rep(1,k))

#parameters to follow for output posterior summaries
params.to.monitor <- c("eta","mu","tau","clust")

#run jags for the bivariate four class model.
jags.fit.CF2x4 <- jags(data=data.package.CF2x4,inits=init.generator,parameters.to.save=params.to.monitor,n.iter=20000,model.file="JAGS_bivar4.txt",n.chains = 4,n.burnin = 4000,n.thin=5 )

#Extract MCMC chains
jagsfitCF2x4.mcmc <- as.mcmc(jags.fit.CF2x4) 

summary(jagsfitCF2x4.mcmc)

#plot trace plots for each recorded parameter to ensure stability of chains

plot(jagsfitCF2x4.mcmc[,"eta[1]"])
plot(jagsfitCF2x4.mcmc[,"eta[2]"])
plot(jagsfitCF2x4.mcmc[,"eta[3]"])
plot(jagsfitCF2x4.mcmc[,"eta[4]"])
plot(jagsfitCF2x4.mcmc[,"mu[1,1]"])
plot(jagsfitCF2x4.mcmc[,"mu[1,2]"])
plot(jagsfitCF2x4.mcmc[,"mu[2,1]"])
plot(jagsfitCF2x4.mcmc[,"mu[2,2]"])
plot(jagsfitCF2x4.mcmc[,"mu[3,1]"])
plot(jagsfitCF2x4.mcmc[,"mu[3,2]"])
plot(jagsfitCF2x4.mcmc[,"mu[4,1]"])
plot(jagsfitCF2x4.mcmc[,"mu[4,2]"])
plot(jagsfitCF2x4.mcmc[,"tau[1,1]"])
plot(jagsfitCF2x4.mcmc[,"tau[2,1]"])
plot(jagsfitCF2x4.mcmc[,"tau[1,2]"])
plot(jagsfitCF2x4.mcmc[,"tau[2,2]"])

#densityplot(jagsfitDL2.mcmc,parameters=c("eta","mu","tau"))

DIC_CF2x4 <- jags.fit.CF2x4$BUGSoutput$DIC


#==================================================================##
#DICHOTIC LISTENING (four class, both time points)

data.package.DL2x4 <- list(N=dim(LI_both_DL)[1],D=dim(LI_both_DL)[2], y=LI_both_DL, k=4,
                           S2= diag(D)*1000, S3= diag(D)*1000, mu_0A=rep(-75, D), mu_0B=rep(-25, D), mu_0C=rep(25, D), mu_0D=rep(75, D), alpha = rep(1,k))

#set parameters to monitor
params.to.monitor <- c("eta","mu","tau","clust")

# Run JAGS model
jags.fit.DL2x4 <- jags(data=data.package.DL2x4,inits=init.generator,parameters.to.save=params.to.monitor,n.iter=20000,model.file="JAGS_bivar4.txt",n.chains = 4,n.burnin = 4000,n.thin=5 )

# Extract MCMC chains
jagsfitDL2x4.mcmc <- as.mcmc(jags.fit.DL2x4) 

summary(jagsfitDL2x4.mcmc)

# Plot MCMC trace plots to ensure stability posterior

plot(jagsfitDL2x4.mcmc[,"eta[1]"])
plot(jagsfitDL2x4.mcmc[,"eta[2]"])
plot(jagsfitDL2x4.mcmc[,"eta[3]"])
plot(jagsfitDL2x4.mcmc[,"eta[4]"])
plot(jagsfitDL2x4.mcmc[,"mu[1,1]"])
plot(jagsfitDL2x4.mcmc[,"mu[1,2]"])
plot(jagsfitDL2x4.mcmc[,"mu[2,1]"])
plot(jagsfitDL2x4.mcmc[,"mu[2,2]"])
plot(jagsfitDL2x4.mcmc[,"mu[3,1]"])
plot(jagsfitDL2x4.mcmc[,"mu[3,2]"])
plot(jagsfitDL2x4.mcmc[,"mu[4,1]"])
plot(jagsfitDL2x4.mcmc[,"mu[4,2]"])
plot(jagsfitDL2x4.mcmc[,"tau[1,1]"])
plot(jagsfitDL2x4.mcmc[,"tau[2,1]"])
plot(jagsfitDL2x4.mcmc[,"tau[1,2]"])
plot(jagsfitDL2x4.mcmc[,"tau[2,2]"])

#densityplot(jagsfitDL2.mcmc,parameters=c("eta","mu","tau"))

DIC_DL2x4 <- jags.fit.DL2x4$BUGSoutput$DIC

#==================================================================##
# Finger tapping ( four class, two time points)
#==================================================================##

#set up the data frame for finger tapping
data.package.FT2x4 <- list(N=dim(LI_both_FT)[1],D=dim(LI_both_FT)[2], y=LI_both_FT, k=4,
                           S2= diag(D)*1000, S3= diag(D)*1000, mu_0A=rep(-75, D), mu_0B=rep(25, D), mu_0C=rep(25, D), mu_0D=rep(75, D), alpha = rep(1,k))

#set parameters to monitor
params.to.monitor <- c("eta","mu","tau","clust")

# Run JAGS model for the finger tapping data.
jags.fit.FT2x4 <- jags(data=data.package.FT2x4,inits=init.generator,parameters.to.save=params.to.monitor,n.iter=20000,model.file="JAGS_bivar4.txt",n.chains = 4,n.burnin = 4000,n.thin=5 )

# Extract the MCMC chains to check the posterior summaries.
jagsfitFT2x4.mcmc <- as.mcmc(jags.fit.DL2x4) 

summary(jagsfitFT2x4.mcmc)

# plot the MCMC chains as trace plots 

plot(jagsfitFT2x4.mcmc[,"eta[1]"])
plot(jagsfitFT2x4.mcmc[,"eta[2]"])
plot(jagsfitFT2x4.mcmc[,"eta[3]"])
plot(jagsfitFT2x4.mcmc[,"eta[4]"])
plot(jagsfitFT2x4.mcmc[,"mu[1,1]"])
plot(jagsfitFT2x4.mcmc[,"mu[1,2]"])
plot(jagsfitFT2x4.mcmc[,"mu[2,1]"])
plot(jagsfitFT2x4.mcmc[,"mu[2,2]"])
plot(jagsfitFT2x4.mcmc[,"mu[3,1]"])
plot(jagsfitFT2x4.mcmc[,"mu[3,2]"])
plot(jagsfitFT2x4.mcmc[,"mu[4,1]"])
plot(jagsfitFT2x4.mcmc[,"mu[4,2]"])
plot(jagsfitFT2x4.mcmc[,"tau[1,1]"])
plot(jagsfitFT2x4.mcmc[,"tau[2,1]"])
plot(jagsfitFT2x4.mcmc[,"tau[1,2]"])
plot(jagsfitFT2x4.mcmc[,"tau[2,2]"])

#densityplot(jagsfitDL2.mcmc,parameters=c("eta","mu","tau"))

# Deviance information criterion for model comparison.
DIC_FT2x4 <- jags.fit.FT2x4$BUGSoutput$DIC

#==================================================================##

#==================================================================##
# Repeat analysis using pivmet package. This ensures label switching 
# is fixed. The analysis is virtually identical model but the previous
# model above uses starting values and more specific priors.

#==================================================================##
# Chimeric faces
#==================================================================##

# run model with label switching adjustment
res_CF_both4<-piv_MCMC(y = as.matrix(LI_both_CF), k = 4, nMC = 10000, software = "rjags")
rel_CF_both4 <- piv_rel(mcmc=res_CF_both4) #extract relative estimates

piv_plot(y = as.matrix(LI_both_CF), mcmc = res_CF_both4, rel_est = rel_CF_both4, par = "mean", type = "chains")
piv_plot(y = as.matrix(LI_both_CF), mcmc = res_CF_both4, rel_est = rel_CF_both4, par = "weight", type = "chains")
#> Description: traceplot of the raw MCMC chains and the relabelled chains for the means parameters (coordinate 1 and 2). Each colored chain corresponds to one of the k distinct parameters of the mixture model. Overlapping chains may reveal that the MCMC sample is not able to distinguish between the components.
piv_plot(y = as.matrix(LI_both_CF), mcmc = res_CF_both4, rel_est = rel_CF_both4, type = "hist")
#> Description: 3d histogram of the data along with the posterior estimates of the relabelled means (triangle points)

#==================================================================##
#Dichotic Listening
#==================================================================##
# run model with label switching adjustment
res_DL_both4<-piv_MCMC(y = as.matrix(LI_both_DL), k = 4, nMC = 10000, software = "rjags")
rel_DL_both4 <- piv_rel(mcmc=res_DL_both4)#extract relative estimates

piv_plot(y = as.matrix(LI_both_DL), mcmc = res_DL_both4, rel_est = rel_DL_both4, par = "mean", type = "chains")
piv_plot(y = as.matrix(LI_both_DL), mcmc = res_DL_both4, rel_est = rel_DL_both4, par = "weight", type = "chains")
#> Description: traceplot of the raw MCMC chains and the relabelled chains for the means parameters (coordinate 1 and 2). Each colored chain corresponds to one of the k distinct parameters of the mixture model. Overlapping chains may reveal that the MCMC sample is not able to distinguish between the components.
piv_plot(y = as.matrix(LI_both_DL), mcmc = res_DL_both4, rel_est = rel_DL_both4, type = "hist")
#> Description: 3d histogram of the data along with the posterior estimates of the relabelled means (triangle points)

#==================================================================##
#Finger tapping
#==================================================================##
# run model with label switching adjustment
res_FT_both4<-piv_MCMC(y = as.matrix(LI_both_FT), k = 4, nMC = 10000, software = "rjags")
rel_FT_both4 <- piv_rel(mcmc=res_FT_both2)#extract relative estimates

piv_plot(y = as.matrix(LI_both_FT), mcmc = res_FT_both4, rel_est = rel_FT_both4, par = "mean", type = "chains")
piv_plot(y = as.matrix(LI_both_FT), mcmc = res_FT_both4, rel_est = rel_FT_both4, par = "weight", type = "chains")
#> Description: traceplot of the raw MCMC chains and the relabelled chains for the means parameters (coordinate 1 and 2). Each colored chain corresponds to one of the k distinct parameters of the mixture model. Overlapping chains may reveal that the MCMC sample is not able to distinguish between the components.
piv_plot(y = as.matrix(LI_both_FT), mcmc = res_FT_both4, rel_est = rel_FT_both4, type = "hist")
#> Description: 3d histogram of the data along with the posterior estimates of the relabelled means (triangle points)

#==================================================================##
#==================================================================##
#==================================================================##
#Model comparison and optimal choice.

model.choice<-c("2-class","3-class","4-class")

CF_DICs <- c(DIC_CF2x2,DIC_CF2x3,DIC_CF2x4)
DL_DICs <- c(DIC_DL2x2,DIC_DL2x3,DIC_DL2x4)
FT_DICs <- c(DIC_FT2x2,DIC_FT2x3,DIC_FT2x4)

compare_DIC_CF <- paste0(model.choice[which.min(CF_DICs)],"_CF")
compare_DIC_DL <- paste0(model.choice[which.min(DL_DICs)],"_DL")
compare_DIC_FT <- paste0(model.choice[which.min(FT_DICs)],"_FT")