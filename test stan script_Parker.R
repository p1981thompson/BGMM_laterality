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


#==================================================================##
# Stan model
#==================================================================##

mix_univ <-"
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
                  iter = 50000)

#==================================================================##

Parker_univ2_DL <- print(fit_univ2_DL, pars = c("mu","sigma","eta"))

traceplot(fit_univ2_DL)
#==================================================================##

library(loo)
library(flextable)
log_lik1_DL<-extract_log_lik(fit_univ2_DL, merge_chains=FALSE)

waic1_DL<-waic(log_lik1_DL)

waic1_DL_tab <- flextable(waic1)

r_eff_DL <- relative_eff(exp(log_lik1_DL), cores = 2)
loo_1_DL <- loo(log_lik1_DL, r_eff = r_eff_DL, cores = 2)

loo_1_DL_tab<-flextable(loo_1_DL)

#==================================================================##
 k <- 2
 nMC <- 50000
 
 res2_DL <- piv_MCMC(y = myLI, k = k, nMC = 10000, 
                  software = "rstan")
 rel2_DL <- piv_rel(res2_DL)
 
 png(filename = "Parker_univ2_DL_trace.png", width = 10, height = 4, units = "in")
 piv_plot(y=myLI, res2_DL, rel2_DL, par = c("mean"), type="chains")
 dev.off()
 
 posterior2_DL <- as.array(res2_DL$stanfit)
 
 png(filename = "Parker_univ2_DL_forest.png", width = 10, height = 4, units = "in")
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
                   iter = 50000)

#==================================================================##

Parker_univ3_DL <- print(fit_univ3_DL, pars = c("mu","sigma","eta"))

traceplot(fit_univ3_DL)

#==================================================================##


log_lik2_DL<-extract_log_lik(fit_univ3_DL, merge_chains=FALSE)

waic2<-waic(log_lik2_DL)

waic2_DL_tab<-flextable(waic2)

r_eff2 <- relative_eff(exp(log_lik2_DL), cores = 2)
loo_2 <- loo(log_lik2_DL, r_eff = r_eff2_DL, cores = 2)

loo_2_DL_tab<-flextable(loo_2_DL)

# Compare
comp1 <- loo_compare(loo_1_DL, loo_2_DL)
print(comp1)

#==================================================================##
k <- 3
nMC <- 50000

res3_DL <- piv_MCMC(y = myLI, k = k, nMC = 10000, 
                 software = "rstan")
rel3_DL <- piv_rel(res3_DL)

png(filename = "Parker_univ3_DL_trace.png", width = 10, height = 4, units = "in")
piv_plot(y=y, res3_DL, rel3_DL, par = c("mean"), type="chains")
dev.off()

posterior3_DL <- as.array(res3_DL$stanfit)

png(filename = "Parker_univ3_DL_forest.png", width = 10, height = 4, units = "in")
mcmc_intervals(posterior3_DL, regex_pars = c("mu","sigma"))
dev.off()

#==================================================================##
# 4 component Gaussian mixture model via Stan
#==================================================================##

dataC = list(N=length(myLI), y=myLI, k=4,
             mu_0=0, B0inv=0.1,
             mu_sigma=0, tau_sigma=2)

#==================================================================##

fit_univ4_DL <-  stan(model_code = mix_univ,
                   data=dataC,
                   chains = 4,
                   iter = 50000)

#==================================================================##

Parker_univ4_DL <- print(fit_univ4_DL, pars = c("mu","sigma","eta"))

traceplot(fit_univ4_DL)

#==================================================================##

library(loo)
log_lik3_DL<-extract_log_lik(fit_univ4_DL, merge_chains=FALSE)

waic3<-waic(log_lik3_DL)
wiac3_DL_tab<-flextable(waic3)

r_eff3_DL <- relative_eff(exp(log_lik3_DL), cores = 2)
loo_3_DL <- loo(log_lik3_DL, r_eff = r_eff3_DL, cores = 2)

loo_3_DL_tab<-flextable(loo_3_DL)

# Compare
comp2 <- loo_compare(loo_1_DL, loo_2_DL, loo_3_DL)
print(comp2)

#==================================================================##
k <- 4
nMC <- 50000

res4_DL <- piv_MCMC(y = myLI, k = k, nMC = 10000, 
                 software = "rstan")
rel4_DL <- piv_rel(res4_DL)

png(filename = "Parker_univ4_DL_trace.png", width = 10, height = 4, units = "in")
piv_plot(y=y, res4_DL, rel4_DL, par = c("mean"), type="chains")
dev.off()

posterior4_DL <- as.array(res4_DL$stanfit)

png(filename = "Parker_univ4_DL_trace.png", width = 10, height = 4, units = "in")
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
                      iter = 50000)

#==================================================================##

Parker_univ2_CF <- print(fit_univ2_CF, pars = c("mu","sigma","eta"))

traceplot(fit_univ2_CF)
#==================================================================##

library(loo)
log_lik1_CF<-extract_log_lik(fit_univ2_CF, merge_chains=FALSE)

waic1<-waic(log_lik1_CF)

waic1_CF_tab<-flextable(waic1)

r_eff_CF <- relative_eff(exp(log_lik1_CF), cores = 2)
loo_1_CF <- loo(log_lik1_CF, r_eff = r_eff_CF, cores = 2)

loo_1_CF_tab<-flextable(loo_1_CF)

print(loo_1_CF)

#==================================================================##
k <- 2
nMC <- 50000

res2_CF <- piv_MCMC(y = myLI2, k = k, nMC = 10000, 
                    software = "rstan")
rel2_CF <- piv_rel(res2_CF)

png(filename = "Parker_univ2_CF_trace.png", width = 10, height = 4, units = "in")
piv_plot(y=myLI, res2_CF, rel2_CF, par = c("mean"), type="chains")
dev.off()

posterior2_CF <- as.array(res2_CF$stanfit)

png(filename = "Parker_univ2_CF_forest.png", width = 10, height = 4, units = "in")
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
                      iter = 50000)

#==================================================================##

Parker_univ3_CF <- print(fit_univ3_CF, pars = c("mu","sigma","eta"))

traceplot(fit_univ3_CF)

#==================================================================##


log_lik2_CF<-extract_log_lik(fit_univ3_CF, merge_chains=FALSE)

waic2<-waic(log_lik2_CF)

waic2_CF_tab<-flextable(waic2)

r_eff2 <- relative_eff(exp(log_lik2_CF), cores = 2)
loo_2 <- loo(log_lik2_CF, r_eff = r_eff2_CF, cores = 2)

loo_2_CF<-flextable(loo_2)



# Compare
comp1 <- loo_compare(loo_1_CF, loo_2_CF)
print(comp1)

#==================================================================##
k <- 3
nMC <- 50000

res3_CF <- piv_MCMC(y = myLI2, k = k, nMC = 10000, 
                    software = "rstan")
rel3_CF <- piv_rel(res3_CF)

png(filename = "Parker_univ3_CF_trace.png", width = 10, height = 4, units = "in")
piv_plot(y=y, res3_CF, rel3_CF, par = c("mean"), type="chains")
dev.off()

posterior3_CF <- as.array(res3_CF$stanfit)

png(filename = "Parker_univ3_CF_forest.png", width = 10, height = 4, units = "in")
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
                      iter = 50000)

#==================================================================##

Parker_univ4_CF <- print(fit_univ4_CF, pars = c("mu","sigma","eta"))

traceplot(fit_univ4_CF)

#==================================================================##

library(loo)
log_lik3_CF<-extract_log_lik(fit_univ4_CF, merge_chains=FALSE)

waic3<-waic(log_lik3_CF)

waic3_CF_tab<-flextable(waic3)

r_eff3_CF <- relative_eff(exp(log_lik3_CF), cores = 2)
loo_3_CF <- loo(log_lik3_CF, r_eff = r_eff3_CF, cores = 2)

loo_3_CF<-flextable(loo_3_CF)


# Compare
comp2 <- loo_compare(loo_1_CF, loo_2_CF, loo_3_CF)
print(comp2)

#==================================================================##
k <- 4
nMC <- 50000

res4_CF <- piv_MCMC(y = myLI2, k = k, nMC = 10000, 
                    software = "rstan")
rel4_CF <- piv_rel(res4_CF)

png(filename = "Parker_univ4_CF_trace.png", width = 10, height = 4, units = "in")
piv_plot(y=y, res4_CF, rel4_CF, par = c("mean"), type="chains")
dev.off()

posterior4_CF <- as.array(res4_CF$stanfit)

png(filename = "Parker_univ4_CF_forest.png", width = 10, height = 4, units = "in")
mcmc_intervals(posterior4_CF, regex_pars = c("mu","sigma"))
dev.off()