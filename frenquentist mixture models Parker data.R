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

#==================================================================##

# Bootstrap likelihood ratio to test which no. of components is best.

LI.DL.boot <- boot.comp(myLI,max.comp=10,mix.type="normalmix", maxit=400,epsilon=1e-2)

LI.CF.boot <- boot.comp(myLI2,max.comp=10,mix.type="normalmix", maxit=400,epsilon=1e-2)

LI.FT.boot <- boot.comp(myLI3,max.comp=10,mix.type="normalmix", maxit=400,epsilon=1e-2)

#==================================================================##
#Checking calibration of mixture

pnormmix <- function(x,mixture) {
  lambda <- mixture$lambda
  k <- length(lambda)
  pnorm.from.mix <- function(x,component) {
    lambda[component]*pnorm(x,mean=mixture$mu[component],
                            sd=mixture$sigma[component])
  }
  pnorms <- sapply(1:k,pnorm.from.mix,x=x)
  return(rowSums(pnorms))
}

#Dichotic Listening
distinct.DL <- sort(unique(myLI))
tcdfs <- pnormmix(distinct.DL,mixture=LI_fit_DL)
ecdfs <- ecdf(myLI)(distinct.DL)
plot(tcdfs,ecdfs,xlab="Theoretical CDF",ylab="Empirical CDF",xlim=c(0,1),
     ylim=c(0,1))
abline(0,1)

#Chimeric faces
distinct.CF <- sort(unique(myLI2))
tcdfs <- pnormmix(distinct.CF,mixture=LI_fit_CF)
ecdfs <- ecdf(myLI2)(distinct.CF)
plot(tcdfs,ecdfs,xlab="Theoretical CDF",ylab="Empirical CDF",xlim=c(0,1),
     ylim=c(0,1))
abline(0,1)

#Finger Tapping
distinct.FT <- sort(unique(myLI3))
tcdfs <- pnormmix(distinct.FT,mixture=LI_fit_FT)
ecdfs <- ecdf(myLI3)(distinct.FT)
plot(tcdfs,ecdfs,xlab="Theoretical CDF",ylab="Empirical CDF",xlim=c(0,1),
     ylim=c(0,1))
abline(0,1)

#==================================================================##
# Sequentially fit models for k=2,3,4

#Dichotic Listening#

LI.DL.k2 <- normalmixEM(myLI,k=2,maxit=100,epsilon=0.01)
LI.DL.k3 <- normalmixEM(myLI,k=3,maxit=100,epsilon=0.01)
LI.DL.k4 <- normalmixEM(myLI,k=4,maxit=100,epsilon=0.01)

#plot density estimates for DL (k=2)
plot_mix_comps <- function(x, mu, sigma, lam) {
  lam * dnorm(x, mu, sigma)
}
data.frame(x = LI.DL.k2$x) %>%
  ggplot() +
  geom_histogram(aes(x, ..density..), binwidth = 2, colour = "grey", 
                 fill = "white") +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(LI.DL.k2$mu[1], LI.DL.k2$sigma[1], lam = LI.DL.k2$lambda[1]),
                colour = "red", lwd = 1.5) +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(LI.DL.k2$mu[2], LI.DL.k2$sigma[2], lam = LI.DL.k2$lambda[2]),
                colour = "blue", lwd = 1.5) +
  ylab("Density")+theme_bw()


#plot density estimates for DL (k=3)
data.frame(x = LI.DL.k3$x) %>%
  ggplot() +
  geom_histogram(aes(x, ..density..), binwidth = 2, colour = "grey", 
                 fill = "white") +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(LI.DL.k3$mu[1], LI.DL.k3$sigma[1], lam = LI.DL.k3$lambda[1]),
                colour = "red", lwd = 1.5) +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(LI.DL.k3$mu[2], LI.DL.k3$sigma[2], lam = LI.DL.k3$lambda[2]),
                colour = "blue", lwd = 1.5) +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(LI.DL.k3$mu[3], LI.DL.k3$sigma[3], lam = LI.DL.k3$lambda[3]),
                colour = "green", lwd = 1.5) +
  ylab("Density")+theme_bw()


#plot density estimates for DL (k=4)
data.frame(x = LI.DL.k4$x) %>%
  ggplot() +
  geom_histogram(aes(x, ..density..), binwidth = 2, colour = "grey", 
                 fill = "white") +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(LI.DL.k4$mu[1], LI.DL.k4$sigma[1], lam = LI.DL.k4$lambda[1]),
                colour = "red", lwd = 1.5) +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(LI.DL.k4$mu[2], LI.DL.k4$sigma[2], lam = LI.DL.k4$lambda[2]),
                colour = "blue", lwd = 1.5) +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(LI.DL.k4$mu[3], LI.DL.k4$sigma[3], lam = LI.DL.k4$lambda[3]),
                colour = "green", lwd = 1.5) +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(LI.DL.k4$mu[4], LI.DL.k4$sigma[4], lam = LI.DL.k4$lambda[4]),
                colour = "purple", lwd = 1.5) +
  ylab("Density")+theme_bw()

#==================================================================##
# Sequentially fit models for k=2,3,4


#Chimeric faces#

LI.CF.k2 <- normalmixEM(myLI2,k=2,maxit=100,epsilon=0.01)
LI.CF.k3 <- normalmixEM(myLI2,k=3,maxit=100,epsilon=0.01)
LI.CF.k4 <- normalmixEM(myLI2,k=4,maxit=100,epsilon=0.01)

#plot density estimates for CF (k=2)
data.frame(x = LI.CF.k2$x) %>%
  ggplot() +
  geom_histogram(aes(x, ..density..), binwidth = 2, colour = "grey", 
                 fill = "white") +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(LI.CF.k2$mu[1], LI.CF.k2$sigma[1], lam = LI.CF.k2$lambda[1]),
                colour = "red", lwd = 1.5) +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(LI.CF.k2$mu[2], LI.CF.k2$sigma[2], lam = LI.CF.k2$lambda[2]),
                colour = "blue", lwd = 1.5) +
  ylab("Density")+theme_bw()

#plot density estimates for CF (k=3)
data.frame(x = LI.CF.k3$x) %>%
  ggplot() +
  geom_histogram(aes(x, ..density..), binwidth = 2, colour = "grey", 
                 fill = "white") +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(LI.CF.k3$mu[1], LI.CF.k3$sigma[1], lam = LI.CF.k3$lambda[1]),
                colour = "red", lwd = 1.5) +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(LI.CF.k3$mu[2], LI.CF.k3$sigma[2], lam = LI.CF.k3$lambda[2]),
                colour = "blue", lwd = 1.5) +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(LI.CF.k3$mu[3], LI.CF.k3$sigma[3], lam = LI.CF.k3$lambda[3]),
                colour = "green", lwd = 1.5) +
  ylab("Density")+theme_bw()

#plot density estimates for CF (k=4)
data.frame(x = LI.CF.k4$x) %>%
  ggplot() +
  geom_histogram(aes(x, ..density..), binwidth = 2, colour = "grey", 
                 fill = "white") +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(LI.CF.k4$mu[1], LI.CF.k4$sigma[1], lam = LI.CF.k4$lambda[1]),
                colour = "red", lwd = 1.5) +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(LI.CF.k4$mu[2], LI.CF.k4$sigma[2], lam = LI.CF.k4$lambda[2]),
                colour = "blue", lwd = 1.5) +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(LI.CF.k4$mu[3], LI.CF.k4$sigma[3], lam = LI.CF.k4$lambda[3]),
                colour = "green", lwd = 1.5) +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(LI.CF.k4$mu[4], LI.CF.k4$sigma[4], lam = LI.CF.k4$lambda[4]),
                colour = "purple", lwd = 1.5) +
  ylab("Density")+theme_bw()

#==================================================================##
# Sequentially fit models for k=2,3,4


#Finger Tapping#

LI.FT.k2 <- normalmixEM(myLI3,k=2,maxit=100,epsilon=0.01)
LI.FT.k3 <- normalmixEM(myLI3,k=3,maxit=100,epsilon=0.01)
LI.FT.k4 <- normalmixEM(myLI3,k=4,maxit=100,epsilon=0.01)

#plot density estimates for FT (k=2)
data.frame(x = LI.FT.k2$x) %>%
  ggplot() +
  geom_histogram(aes(x, ..density..), binwidth = 2, colour = "grey", 
                 fill = "white") +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(LI.FT.k2$mu[1], LI.FT.k2$sigma[1], lam = LI.FT.k2$lambda[1]),
                colour = "red", lwd = 1.5) +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(LI.FT.k2$mu[2], LI.FT.k2$sigma[2], lam = LI.FT.k2$lambda[2]),
                colour = "blue", lwd = 1.5) +
  ylab("Density")+theme_bw()

#plot density estimates for FT (k=3)
data.frame(x = LI.FT.k3$x) %>%
  ggplot() +
  geom_histogram(aes(x, ..density..), binwidth = 2, colour = "grey", 
                 fill = "white") +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(LI.FT.k3$mu[1], LI.FT.k3$sigma[1], lam = LI.FT.k3$lambda[1]),
                colour = "red", lwd = 1.5) +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(LI.FT.k3$mu[2], LI.FT.k3$sigma[2], lam = LI.FT.k3$lambda[2]),
                colour = "blue", lwd = 1.5) +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(LI.CF.k3$mu[3], LI.CF.k3$sigma[3], lam = LI.CF.k3$lambda[3]),
                colour = "green", lwd = 1.5) +
  ylab("Density")+theme_bw()

#plot density estimates for FT (k=4)
data.frame(x = LI.FT.k4$x) %>%
  ggplot() +
  geom_histogram(aes(x, ..density..), binwidth = 2, colour = "grey", 
                 fill = "white") +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(LI.FT.k4$mu[1], LI.FT.k4$sigma[1], lam = LI.FT.k4$lambda[1]),
                colour = "red", lwd = 1.5) +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(LI.FT.k4$mu[2], LI.FT.k4$sigma[2], lam = LI.FT.k4$lambda[2]),
                colour = "blue", lwd = 1.5) +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(LI.FT.k4$mu[3], LI.FT.k4$sigma[3], lam = LI.FT.k4$lambda[3]),
                colour = "green", lwd = 1.5) +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(LI.FT.k4$mu[4], LI.FT.k4$sigma[4], lam = LI.FT.k4$lambda[4]),
                colour = "purple", lwd = 1.5) +
  ylab("Density")+theme_bw()


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

# Run for each variable to see which k component model is selected.

CVfun(myLI)
CVfun(myLI2)
CVfun(myLI3)

#Pretty inconclusive! mostly pointing toward 2 group mixtures

#======================================================================#

# Now test whether we have homogeneity in the data, i.e. whether this is a mixture at all or just a single group.

MixtureInf::emtest.norm(myLI, niter = 500) #DL looks like 2 groups
MixtureInf::emtest.norm(myLI2, niter = 500) #DL looks like 2 groups
MixtureInf::emtest.norm(myLI3, niter = 500) #DL looks like 2 groups

#======================================================================================#
# Bivariate mixture models using both time points.
#======================================================================================#

mu <- list(c(-50, 50), c(-50, 50)) #starting values for the EM algorithm
DL_bivarEM<-mvnormalmixEM(LI_both_DL, k = 2, mu=mu, arbmean = TRUE,epsilon = 1e-02)
plot(DL_bivarEM, density = TRUE,  marginal = TRUE)

boot.comp(y=as.matrix(LI_both_DL),max.comp=6,arbmean=TRUE, mix.type = "mvnormalmix", hist=TRUE)



# I'm not really convinced by the following for CF. I personally think this is one group.

mu <- list(c(-20, 20), c(-20, 20))#starting values for the EM algorithm
CF_bivarEM<-mvnormalmixEM(LI_both_CF, k = 2, mu=mu, arbvar = TRUE,epsilon = 1e-02)
plot(CF_bivarEM, density = TRUE,  marginal = TRUE)

boot.comp(y=as.matrix(LI_both_CF),max.comp=6,arbmean=TRUE, mix.type = "mvnormalmix", hist=TRUE)

#This appears to have the same location but different scales, so changed to "arbvar" argument.

mu <- list(c(-20, 20), c(-20, 20))#starting values for the EM algorithm
FT_bivarEM<-mvnormalmixEM(LI_both_FT, k = 2, mu=mu, arbvar = TRUE,epsilon = 1e-02)
plot(FT_bivarEM, density = TRUE,  marginal = TRUE)

boot.comp(y=as.matrix(LI_both_FT),max.comp=6,arbmean=TRUE, mix.type = "mvnormalmix", hist=TRUE)