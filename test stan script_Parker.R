library(pivmet)
library(bayesmix)
library(bayesplot)
library(tidyverse)
library(rstan)


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
#
# Stan model
#==================================================================##

mix_univ_new_priors <-"
        data {
          int<lower=1> k;          // number of mixture components
          int<lower=1> N;          // number of data points
          real y[N];               // observations
          real mu_0;               // mean hyperparameter
          real<lower=0> B0inv;     // mean hyperprecision
          real mu_sigma;           // sigma hypermean
          real<lower=0> tau_sigma; // sigma hyper sd
          }
        parameters {
          simplex[k] eta;             // mixing proportions
          ordered[k] mu;              // locations of mixture components - means don't switch.
          vector<lower=0>[k] sigma;   // scales of mixture components
        }
          
      transformed parameters{
          vector[k] log_eta = log(eta);  // cache log calculation
          vector[k] pz[N];
          simplex[k] exp_pz[N];
              for (n in 1:N){
                  pz[n] =   normal_lpdf(y[n]|mu, sigma)+
                            log_eta-
                            log_sum_exp(normal_lpdf(y[n]|mu, sigma)+
                            log_eta);
                  exp_pz[n] = exp(pz[n]);
                            }
          }
      model {
        sigma ~ lognormal(mu_sigma, tau_sigma);
        mu ~ normal(mu_0, 1/B0inv);
            for (n in 1:N) {
              vector[k] lps = log_eta;
                for (j in 1:k){
                    lps[j] += normal_lpdf(y[n] | mu[j], sigma[j]);
                    target+=pz[n,j];
                    }
              target += log_sum_exp(lps);
                  }
          }
     generated quantities{
        int<lower=1, upper=k> z[N];
        vector[N] log_lik;
        vector[k] lps;
          for (n in 1:N){
              z[n] = categorical_rng(exp_pz[n]);
              
                for (j in 1:k){
                    lps[j] = normal_lpdf(y[n] | mu[j], sigma[j]);
                    }
              log_lik[n] = log_sum_exp(lps);
            }
      }
      "
#==================================================================##
# Dichotic listening
#==================================================================##
# 2 component Gaussian mixture model via Stan
#==================================================================##
#
data = list(N=length(myLI), y=myLI, k=2,
            mu_0=0, B0inv=0.1,
            mu_sigma=0, tau_sigma=2)

#==================================================================##

fit_univ2_DL <-  stan(model_code = mix_univ,
                      data=data,
                      chains = 4,
                      iter = 50000,
                      thin=5)

#==================================================================##

Parker_univ2_DL <- print(fit_univ2_DL, pars = c("mu","sigma","eta"))

png(filename = "Parker_univ2_DL_trace_stan_direct.png", width = 12, height = 6, units = "in",res=300)
traceplot(fit_univ2_DL)
dev.off()
#==================================================================##

library(loo)
library(flextable)
log_lik1_DL<-extract_log_lik(fit_univ2_DL, merge_chains=FALSE)

waic1_DL<-waic(log_lik1_DL)

r_eff_DL <- relative_eff(exp(log_lik1_DL), cores = 2)
loo_1_DL <- loo(log_lik1_DL, r_eff = r_eff_DL, cores = 2)


#==================================================================##
k <- 2
nMC <- 10000

res2_DL <- piv_MCMC(y = myLI, k = k, nMC = nMC, 
                    software = "rstan")
rel2_DL <- piv_rel(res2_DL)

png(filename = "Parker_univ2_DL_trace.png", width = 10, height = 4, units = "in",res=300)
piv_plot(y=myLI, res2_DL, rel2_DL, par = c("mean"), type="chains")
dev.off()

png(filename = "Parker_univ2_DL_hist.png", width = 10, height = 4, units = "in",res=300)
piv_plot(y=myLI, res2_DL, rel2_DL, par = c("mean"), type="hist")
dev.off()

posterior2_DL <- as.array(res2_DL$stanfit)

png(filename = "Parker_univ2_DL_forest.png", width = 10, height = 4, units = "in",res=300)
mcmc_intervals(posterior2_DL, regex_pars = c("mu","sigma"))
dev.off()
#==================================================================##
# 3 component Gaussian mixture model via Stan
#==================================================================##

dataB = list(N=length(myLI), y=myLI, k=3,
             mu_0=0, B0inv=0.1,
             mu_sigma=0, tau_sigma=2)

#==================================================================##

fit_univ3_DL <-  stan(model_code = mix_univ,
                      data=dataB,
                      chains = 4,
                      iter = 50000,
                      thin=5)

#==================================================================##

Parker_univ3_DL <- print(fit_univ3_DL, pars = c("mu","sigma","eta"))

png(filename = "Parker_univ3_DL_trace_stan_direct.png", width = 12, height = 6, units = "in",res=300)
traceplot(fit_univ3_DL, pars = c("mu","sigma","eta"))
dev.off()
#==================================================================##


log_lik2_DL<-extract_log_lik(fit_univ3_DL, merge_chains=FALSE)

waic2<-waic(log_lik2_DL)


r_eff2_DL <- relative_eff(exp(log_lik2_DL), cores = 2)
loo_2_DL <- loo(log_lik2_DL, r_eff = r_eff2_DL, cores = 2)


# Compare
comp1 <- loo_compare(loo_1_DL, loo_2_DL)
print(comp1)

#==================================================================##
k <- 3
nMC <- 10000

res3_DL <- piv_MCMC(y = myLI, k = k, nMC = nMC, 
                    software = "rstan")
rel3_DL <- piv_rel(res3_DL)

png(filename = "Parker_univ3_DL_trace.png", width = 10, height = 4, units = "in",res=300)
piv_plot(y=myLI, res3_DL, rel3_DL, par = c("mean"), type="chains")
dev.off()

png(filename = "Parker_univ3_DL_hist.png", width = 10, height = 4, units = "in",res=300)
piv_plot(y=myLI, res3_DL, rel3_DL, par = c("mean"), type="hist")
dev.off()

posterior3_DL <- as.array(res3_DL$stanfit)

png(filename = "Parker_univ3_DL_forest.png", width = 10, height = 4, units = "in",res=300)
mcmc_intervals(posterior3_DL, regex_pars = c("mu","sigma"))
dev.off()

#==================================================================##
# 4 component Gaussian mixture model via Stan
#==================================================================##

dataC = list(N=length(myLI), y=myLI, k=4,
             mu_0=0, B0inv=0.02,
             mu_sigma=0, tau_sigma=2)

#==================================================================##

fit_univ4_DL <-  stan(model_code = mix_univ,
                      data=dataC,
                      chains = 4,
                      iter = 50000,
                      thin=5)

#==================================================================##

Parker_univ4_DL <- print(fit_univ4_DL, pars = c("mu","sigma","eta"))

png(filename = "Parker_univ4_DL_trace_stan_direct.png", width = 12, height = 6, units = "in",res=300)
traceplot(fit_univ4_DL, pars = c("mu","sigma","eta"))
dev.off()

#==================================================================##

library(loo)
log_lik3_DL<-extract_log_lik(fit_univ4_DL, merge_chains=FALSE)

waic3<-waic(log_lik3_DL)


r_eff3_DL <- relative_eff(exp(log_lik3_DL), cores = 2)
loo_3_DL <- loo(log_lik3_DL, r_eff = r_eff3_DL, cores = 2)


# Compare
comp2 <- loo_compare(loo_1_DL, loo_2_DL, loo_3_DL)
print(comp2)

#==================================================================##
k <- 4
nMC <- 20000

res4_DL <- piv_MCMC(y = myLI, k = k, nMC = nMC, 
                    software = "rstan",priors = list(mu_sigma = 0, tau_sigma = 3, B0inv=0.01))
rel4_DL <- piv_rel(res4_DL)

png(filename = "Parker_univ4_DL_trace.png", width = 10, height = 4, units = "in",res=300)
piv_plot(y=myLI, res4_DL, rel4_DL, par = c("mean"), type="chains")
dev.off()

png(filename = "Parker_univ4_DL_hist.png", width = 10, height = 4, units = "in",res=300)
piv_plot(y=myLI, res4_DL, rel4_DL, par = c("mean"), type="hist")
dev.off()

posterior4_DL <- as.array(res4_DL$stanfit)

png(filename = "Parker_univ4_DL_forest.png", width = 10, height = 4, units = "in",res=300)
mcmc_intervals(posterior4_DL, regex_pars = c("mu","sigma"))
dev.off()
#==================================================================##
#==================================================================##
# Chimeric faces
#==================================================================##
# 2 component Gaussian mixture model via Stan
#==================================================================##

data = list(N=length(myLI2), y=myLI2, k=2,
            mu_0=0, B0inv=0.1,
            mu_sigma=0, tau_sigma=2)

#==================================================================##

fit_univ2_CF <-  stan(model_code = mix_univ,
                      data=data,
                      chains = 4,
                      iter = 50000,
                      thin=5)

#==================================================================##

Parker_univ2_CF <- print(fit_univ2_CF, pars = c("mu","sigma","eta"))

png(filename = "Parker_univ2_CF_trace_stan_direct.png", width = 12, height = 6, units = "in",res=300)
traceplot(fit_univ2_CF)
dev.off()

#==================================================================##

library(loo)
log_lik1_CF<-extract_log_lik(fit_univ2_CF, merge_chains=FALSE)

waic1_CF<-waic(log_lik1_CF)


r_eff_CF <- relative_eff(exp(log_lik1_CF), cores = 2)
loo_1_CF <- loo(log_lik1_CF, r_eff = r_eff_CF, cores = 2)


#==================================================================##
k <- 2
nMC <- 10000

res2_CF <- piv_MCMC(y = myLI2, k = k, nMC = nMC, 
                    software = "rstan",priors = list(mu_sigma = 0, tau_sigma = 3, B0inv=0.01))
rel2_CF <- piv_rel(res2_CF)

png(filename = "Parker_univ2_CF_trace.png", width = 10, height = 4, units = "in",res=300)
piv_plot(y=myLI2, res2_CF, rel2_CF, par = c("mean"), type="chains")
dev.off()

png(filename = "Parker_univ2_CF_hist.png", width = 10, height = 4, units = "in",res=300)
piv_plot(y=myLI2, res2_CF, rel2_CF, par = c("mean"), type="hist")
dev.off()

posterior2_CF <- as.array(res2_CF$stanfit)

png(filename = "Parker_univ2_CF_forest.png", width = 10, height = 4, units = "in",res=300)
mcmc_intervals(posterior2_CF, regex_pars = c("mu","sigma"))
dev.off()
#==================================================================##
# 3 component Gaussian mixture model via Stan
#==================================================================##

dataB = list(N=length(myLI2), y=myLI2, k=3,
             mu_0=0, B0inv=0.1,
             mu_sigma=0, tau_sigma=2)

#==================================================================##

fit_univ3_CF <-  stan(model_code = mix_univ,
                      data=dataB,
                      chains = 4,
                      iter = 50000,
                      thin=5)

#==================================================================##

Parker_univ3_CF <- print(fit_univ3_CF, pars = c("mu","sigma","eta"))

png(filename = "Parker_univ3_CF_trace_stan_direct.png", width = 12, height = 6, units = "in",res=300)
traceplot(fit_univ3_CF, pars = c("mu","sigma","eta"))
dev.off()
#==================================================================##


log_lik2_CF<-extract_log_lik(fit_univ3_CF, merge_chains=FALSE)

waic2_CF<-waic(log_lik2_CF)


r_eff2_CF <- relative_eff(exp(log_lik2_CF), cores = 2)
loo_2_CF <- loo(log_lik2_CF, r_eff = r_eff2_CF, cores = 2)

# Compare
comp1_CF <- loo_compare(loo_1_CF, loo_2_CF)
print(comp1)

#==================================================================##
k <- 3
nMC <- 10000

res3_CF <- piv_MCMC(y = myLI2, k = k, nMC = nMC, 
                    software = "rstan",priors = list(mu_sigma = 0, tau_sigma = 3, B0inv=0.01))
rel3_CF <- piv_rel(res3_CF)

png(filename = "Parker_univ3_CF_trace.png", width = 10, height = 4, units = "in",res=300)
piv_plot(y=myLI2, res3_CF, rel3_CF, par = c("mean"), type="chains")
dev.off()

png(filename = "Parker_univ3_CF_hist.png", width = 10, height = 4, units = "in",res=300)
piv_plot(y=myLI2, res3_CF, rel3_CF, par = c("mean"), type="hist")
dev.off()

posterior3_CF <- as.array(res3_CF$stanfit)

png(filename = "Parker_univ3_CF_forest.png", width = 10, height = 4, units = "in",res=300)
mcmc_intervals(posterior3_CF, regex_pars = c("mu","sigma"))
dev.off()
#==================================================================##
# 4 component Gaussian mixture model via Stan
#==================================================================##

dataC = list(N=length(myLI2), y=myLI2, k=4,
             mu_0=0, B0inv=0.1,
             mu_sigma=0, tau_sigma=2)

#==================================================================##

fit_univ4_CF <-  stan(model_code = mix_univ,
                      data=dataC,
                      chains = 4,
                      iter = 50000,
                      thin=5)

#==================================================================##

Parker_univ4_CF <- print(fit_univ4_CF, pars = c("mu","sigma","eta"))

png(filename = "Parker_univ4_CF_trace_stan_direct.png", width = 12, height = 6, units = "in",res=300)
traceplot(fit_univ4_CF, pars = c("mu","sigma","eta"))
dev.off()
#==================================================================##

library(loo)
log_lik3_CF<-extract_log_lik(fit_univ4_CF, merge_chains=FALSE)

waic3_CF<-waic(log_lik3_CF)

r_eff3_CF <- relative_eff(exp(log_lik3_CF), cores = 2)
loo_3_CF <- loo(log_lik3_CF, r_eff = r_eff3_CF, cores = 2)

# Compare
comp2_CF <- loo_compare(loo_1_CF, loo_2_CF, loo_3_CF)


#==================================================================##
k <- 4
nMC <- 10000

res4_CF <- piv_MCMC(y = myLI2, k = k, nMC = nMC, 
                    software = "rstan")
write.csv(x=res4_CF$stanfit,file="/Volumes/PSYHOME/PSYRES/pthompson/DVMB/BGMM_laterality/model_est_res4_CF.csv")
rel4_CF <- piv_rel(res4_CF)

png(filename = "Parker_univ4_CF_trace.png", width = 10, height = 4, units = "in",res=300)
piv_plot(y=myLI2, res4_CF, rel4_CF, par = c("mean"), type="chains")
dev.off()

posterior4_CF <- as.array(res4_CF$stanfit)

png(filename = "Parker_univ4_CF_forest.png", width = 10, height = 4, units = "in",res=300)
mcmc_intervals(posterior4_CF, regex_pars = c("mu","sigma"))
dev.off()

save.image("/Volumes/PSYHOME/PSYRES/pthompson/DVMB/BGMM_laterality/feb24.Rdata")
savehistory("/Volumes/PSYHOME/PSYRES/pthompson/DVMB/BGMM_laterality/hist_feb24.Rhistory")