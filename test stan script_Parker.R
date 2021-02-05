library(pivmet)
library(bayesmix)
library(bayesplot)
library(tidyverse)
library(rstan)

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

#==================================================================##
# k <- 3
# nMC <- 15000
# 
# res2 <- piv_MCMC(y = myLI, k = k, nMC = 10000, 
#                  software = "rstan")
# rel2 <- piv_rel(res2)
# piv_plot(y=y, res2, rel2, par = c("mean"), type="chains")
# 
# posterior <- as.array(res2$stanfit)
# mcmc_intervals(posterior, regex_pars = c("mu"))

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
# 2 component Gaussian mixture model via Stan
#==================================================================##

data = list(N=length(myLI), y=myLI, k=2,
            mu_0=0, B0inv=0.1,
            mu_sigma=0, tau_sigma=2)

#==================================================================##

fit_univ2 <-  stan(model_code = mix_univ,
                  data=data,
                  chains = 4,
                  iter = 50000)

#==================================================================##

printed <- cat(print(fit_univ2, pars =c("mu", "eta", "sigma")))
Parker_univ2 <- rstan::extract(fit_univ2)

traceplot(fit_univ2)
#==================================================================##

library(loo)
log_lik1<-extract_log_lik(fit_univ2, merge_chains=FALSE)

waic1<-waic(log_lik1)
r_eff <- relative_eff(exp(log_lik1), cores = 2)
loo_1 <- loo(log_lik1, r_eff = r_eff, cores = 2)
print(loo_1)

#==================================================================##
 k <- 2
 nMC <- 50000
 
 res2 <- piv_MCMC(y = myLI, k = k, nMC = 10000, 
                  software = "rstan")
 rel2 <- piv_rel(res2)
 piv_plot(y=myLI, res2, rel2, par = c("mean"), type="chains")
 
 posterior2 <- as.array(res2$stanfit)
 mcmc_intervals(posterior2, regex_pars = c("mu","sigma"))

#==================================================================##
# 3 component Gaussian mixture model via Stan
#==================================================================##

dataB = list(N=length(myLI), y=myLI, k=3,
             mu_0=0, B0inv=0.1,
             mu_sigma=0, tau_sigma=2)

#==================================================================##

fit_univ3 <-  stan(model_code = mix_univ,
                   data=dataB,
                   chains = 4,
                   iter = 50000)

#==================================================================##

printed <- cat(print(fit_univ3, pars =c("mu", "eta", "sigma")))
Parker_univ3 <- rstan::extract(fit_univ3)

traceplot(fit_univ3)

#==================================================================##

library(loo)
log_lik2<-extract_log_lik(fit_univ3, merge_chains=FALSE)

waic2<-waic(log_lik2)
r_eff2 <- relative_eff(exp(log_lik2), cores = 2)
loo_2 <- loo(log_lik2, r_eff = r_eff2, cores = 2)
print(loo_2)

# Compare
comp <- loo_compare(loo_1, loo_2)
print(comp)

#==================================================================##
k <- 3
nMC <- 50000

res3 <- piv_MCMC(y = myLI, k = k, nMC = 10000, 
                 software = "rstan")
rel3 <- piv_rel(res3)
piv_plot(y=y, res3, rel3, par = c("mean"), type="chains")

posterior3 <- as.array(res3$stanfit)
mcmc_intervals(posterior3, regex_pars = c("mu","sigma"))
#==================================================================##
# 4 component Gaussian mixture model via Stan
#==================================================================##

dataC = list(N=length(myLI), y=myLI, k=4,
             mu_0=0, B0inv=0.1,
             mu_sigma=0, tau_sigma=2)

#==================================================================##

fit_univ4 <-  stan(model_code = mix_univ,
                   data=dataC,
                   chains = 4,
                   iter = 50000)

#==================================================================##

printed <- cat(print(fit_univ4, pars =c("mu", "eta", "sigma")))
Parker_univ4 <- rstan::extract(fit_univ4)

traceplot(fit_univ4)

#==================================================================##

library(loo)
log_lik3<-extract_log_lik(fit_univ4, merge_chains=FALSE)

waic3<-waic(log_lik3)
r_eff3 <- relative_eff(exp(log_lik3), cores = 2)
loo_3 <- loo(log_lik3, r_eff = r_eff3, cores = 2)
print(loo_3)

# Compare
comp <- loo_compare(loo_1, loo_2, loo_3)
print(comp)

#==================================================================##
k <- 4
nMC <- 50000

res4 <- piv_MCMC(y = myLI, k = k, nMC = 10000, 
                 software = "rstan")
rel4 <- piv_rel(res4)
piv_plot(y=y, res4, rel4, par = c("mean"), type="chains")

posterior4 <- as.array(res4$stanfit)
mcmc_intervals(posterior4, regex_pars = c("mu","sigma"))
