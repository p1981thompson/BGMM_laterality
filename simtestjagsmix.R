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

# sim.dat<-rnormmix(1000, lambda=c(0.3,0.7), mu=c(20,-40), sigma=c(10,15))

# sim.dat<-rnormmix(1000, lambda=c(0.3,0.7), mu=c(-10,10), sigma=c(10,10))

sim.dat<-rnormmix(1000, lambda=c(0.3,0.7), mu=c(-10,10), sigma=c(10,10))

init.generator1 <- function(){ list(
  mu = c(runif(1, -40,0),runif(1, 0,40)),
  tau = 1/runif(2, 1,10),
  eta = runif(2, 0,1))
}

data.package.sim2 <- list(N=length(sim.dat), y=sim.dat, k=2,
                         b01=-20,b02=20, B0inv=0.05,
                         nu0Half=1,
                         g0Half=1,g0G0Half=20, e=c(1,1))


params.to.monitor <- c("eta","mu","tau","S0")

jags.fit.sim2 <- jags(data=data.package.sim2,inits=init.generator1,parameters.to.save=params.to.monitor,n.iter=30000,model.file="JAGS_mod1.txt",n.chains = 4,n.burnin = 10000,n.thin=5 )

jagsfitsim2.mcmc <- as.mcmc(jags.fit.sim2) 

summary(jagsfitsim2.mcmc)

plot(jagsfitsim2.mcmc[,"eta[1]"])
plot(jagsfitsim2.mcmc[,"eta[2]"])
plot(jagsfitsim2.mcmc[,"mu[1]"])
plot(jagsfitsim2.mcmc[,"mu[2]"])
plot(jagsfitsim2.mcmc[,"tau[1]"])
plot(jagsfitsim2.mcmc[,"tau[2]"])

DIC_sim2 <- jags.fit.sim2$BUGSoutput$DIC



#================================================================================#


data(fish)
y <- fish[,1]
hist(y, breaks=40, prob = TRUE, cex.lab=1.6,
     main ="Fishery data", cex.main =1.7,
     col="navajowhite1", border="navajowhite1")
lines(density(y), lty=1, lwd=3, col="blue")

k <- 5
nMC <- 15000
res <- piv_MCMC(y = y, k = k, nMC = nMC, 
                burn = 0.5*nMC, software = "rjags")

