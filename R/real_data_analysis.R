# example for the real data application 
#================================================================================================================
# Application of the methods to the hydrological dataset from Istanbul, Türkiye
#
# author: Fatih Kızılaslan (fkizilaslan@yahoo.com) and Gülce Cüran
# date: 22-December-2024
#================================================================================================================

rm(list=ls())
source("functions.R")   
library(Newdistns)
library(MCMCpack)
library(readxl)
library(xtable)
library(goftest)

data                             <- read_excel("real_data.xlsx")
X                                <- data$dam.reserved.water.million.m3 # strength variable
Y                                <- data$consumption.million.m3     # stress variable 

# for data X
seed_number                      <- 12024
set.seed(seed_number)
initials.x                       <- runif(3, 0.1, 2)

g1                               <-  "weibull"
optim.X_BFGS_result.weibull      <- optim(par = initials.x, fn = log_likelihood, gr = grad_log_likelihood,
                                          g = g1, data = X, method = "BFGS" , control = list(maxit = 20000), hessian = T)
BFGS_estimates.X.weibull         <- optim.X_BFGS_result.weibull$par

KS.X.weibull                     <- ks.test(X,  "MO_cdf",   par= BFGS_estimates.X.weibull, g = g1 )
AD.X.weibull                     <- ad.test(X,  "MO_cdf",   par= BFGS_estimates.X.weibull, g = g1 )
CVM.X.weibull                    <- cvm.test(X, "MO_cdf",   par= BFGS_estimates.X.weibull, g = g1 )

test.X.weibull                   <- c(KStest=KS.X.weibull$statistic, KSp=KS.X.weibull$p.value,
                                      ADtest=AD.X.weibull$statistic, ADp=AD.X.weibull$p.value,
                                      CVMtest=CVM.X.weibull$statistic, CVMp=CVM.X.weibull$p.value)
model.X.weibull                  <- Model_evaluation_MOG(g = g1, par=BFGS_estimates.X.weibull, data = X)


g1                               <-  "burrxii"
set.seed(seed_number)
initials.x.burr                  <- runif(3, 0.1, 1)
optim.X_BFGS_result.burrxii      <- optim(par = initials.x.burr, fn = log_likelihood, gr = grad_log_likelihood,
                                     g = g1, data = X, method = "BFGS" , control = list(maxit = 20000), hessian = T)
BFGS_estimates.X.burrxii         <- optim.X_BFGS_result.burrxii$par

KS.X.burrxii                     <- ks.test(X,  "MO_cdf",   par= BFGS_estimates.X.burrxii, g = g1 )
AD.X.burrxii                     <- ad.test(X,  "MO_cdf",   par= BFGS_estimates.X.burrxii, g = g1 )
CVM.X.burrxii                    <- cvm.test(X, "MO_cdf",   par= BFGS_estimates.X.burrxii, g = g1 )
test.X.burrxii                   <- c(KStest=KS.X.burrxii$statistic, KSp=KS.X.burrxii$p.value,
                                      ADtest=AD.X.burrxii$statistic, ADp=AD.X.burrxii$p.value,
                                      CVMtest=CVM.X.burrxii$statistic, CVMp=CVM.X.burrxii$p.value)
model.X.burrxii                  <- Model_evaluation_MOG(g = g1, par=BFGS_estimates.X.burrxii, data = X)


g1                               <-  "chen"
set.seed(seed_number)
initials.x.chen                  <- runif(3, 0.1, 1)
optim.X_BFGS_result.chen         <- optim(par = initials.x.chen, fn = log_likelihood, gr = grad_log_likelihood,
                                       g = g1, data = X, method = "BFGS" , control = list(maxit = 20000), hessian = T)
BFGS_estimates.X.chen            <- optim.X_BFGS_result.chen$par

KS.X.chen                        <- ks.test(X,  "MO_cdf",   par= BFGS_estimates.X.chen, g = g1 )
AD.X.chen                        <- ad.test(X,  "MO_cdf",   par= BFGS_estimates.X.chen, g = g1 )
CVM.X.chen                       <- cvm.test(X, "MO_cdf",   par= BFGS_estimates.X.chen, g = g1 )
test.X.chen                      <- c(KStest=KS.X.chen$statistic, KSp=KS.X.chen$p.value,
                                      ADtest=AD.X.chen$statistic, ADp=AD.X.chen$p.value,
                                      CVMtest=CVM.X.chen$statistic, CVMp=CVM.X.chen$p.value)
model.X.chen                     <- Model_evaluation_MOG(g = g1, par=BFGS_estimates.X.chen, data = X)

data.X.results                   <- rbind(BurrXII = c(BFGS_estimates.X.burrxii, unlist(model.X.burrxii) ),
                                          Chen    = c(BFGS_estimates.X.chen, unlist(model.X.chen) ),
                                          Weibull = c(BFGS_estimates.X.weibull, unlist(model.X.weibull) )
                                          )
colnames( data.X.results)[1:3]   <- c("alpha1", "r1", "s1")
round(data.X.results, 5)
round( rbind(test.X.burrxii, test.X.chen, test.X.weibull), 5)


# for data Y
set.seed(seed_number + 10)
initials.y.weibull               <- runif(3, 0.1, 2)
g2                               <- "weibull"
optim.Y_BFGS_result.weibull      <- optim(par = initials.y.weibull, fn = log_likelihood, gr = grad_log_likelihood,
                                       g = g2, data = Y, method = "BFGS" , control = list(maxit = 20000), hessian = T)
BFGS_estimates.Y.weibull         <- optim.Y_BFGS_result.weibull$par

KS.Y.weibull                     <- ks.test(Y,  "MO_cdf",   par= BFGS_estimates.Y.weibull, g = g2 )
AD.Y.weibull                     <- ad.test(Y,  "MO_cdf",   par= BFGS_estimates.Y.weibull, g = g2 )
CVM.Y.weibull                    <- cvm.test(Y, "MO_cdf",   par= BFGS_estimates.Y.weibull, g = g2 )
test.Y.weibull                   <- c(KStest=KS.Y.weibull$statistic, KSp=KS.Y.weibull$p.value,
                                      ADtest=AD.Y.weibull$statistic, ADp=AD.Y.weibull$p.value,
                                      CVMtest=CVM.Y.weibull$statistic, CVMp=CVM.Y.weibull$p.value)
model.Y.weibull                  <- Model_evaluation_MOG(g = g2, par=BFGS_estimates.Y.weibull, data = Y )


g2                               <- "burrxii"
set.seed(seed_number + 10)
initials.y.burrxii               <- runif(3, 0.1, 1)
optim.Y_BFGS_result.burrxii      <- optim(par = initials.y.burrxii, fn = log_likelihood, gr = grad_log_likelihood,
                                      g = g2, data = Y, method = "BFGS" , control = list(maxit = 20000), hessian = T)
BFGS_estimates.Y.burrxii         <- optim.Y_BFGS_result.burrxii$par

KS.Y.burrxii                     <- ks.test(Y,  "MO_cdf",   par= BFGS_estimates.Y.burrxii, g = g2 )
AD.Y.burrxii                     <- ad.test(Y,  "MO_cdf",   par= BFGS_estimates.Y.burrxii, g = g2 )
CVM.Y.burrxii                    <- cvm.test(Y, "MO_cdf",   par= BFGS_estimates.Y.burrxii, g = g2 )

test.Y.burrxii                   <- c(KStest=KS.Y.burrxii$statistic, KSp=KS.Y.burrxii$p.value,
                                      ADtest=AD.Y.burrxii$statistic, ADp=AD.Y.burrxii$p.value,
                                      CVMtest=CVM.Y.burrxii$statistic, CVMp=CVM.Y.burrxii$p.value)
model.Y.burrxii                  <- Model_evaluation_MOG(g = g2, par=BFGS_estimates.Y.burrxii, data = Y )


g2                               <- "chen"
set.seed(seed_number+10)
initials.y.chen                  <- runif(3,0.1,1)
optim.Y_BFGS_result.chen         <- optim(par = initials.y.chen, fn = log_likelihood, gr = grad_log_likelihood,
                                       g = g2, data = Y, method = "BFGS" , control = list(maxit = 20000), hessian = T)
BFGS_estimates.Y.chen            <- optim.Y_BFGS_result.chen$par

KS.Y.chen                        <- ks.test(Y,  "MO_cdf",   par= BFGS_estimates.Y.chen, g = g2 )
AD.Y.chen                        <- ad.test(Y,  "MO_cdf",   par= BFGS_estimates.Y.chen, g = g2 )
CVM.Y.chen                       <- cvm.test(Y, "MO_cdf",   par= BFGS_estimates.Y.chen, g = g2 )
test.Y.chen                      <- c(KStest=KS.Y.chen$statistic, KSp=KS.Y.chen$p.value,
                                      ADtest=AD.Y.chen$statistic, ADp=AD.Y.chen$p.value,
                                      CVMtest=CVM.Y.chen$statistic, CVMp=CVM.Y.chen$p.value)
model.Y.chen                     <- Model_evaluation_MOG(g = g2, par=BFGS_estimates.Y.chen, data = Y )


data.Y.results                   <- rbind(BurrXII = c(BFGS_estimates.Y.burrxii, unlist(model.Y.burrxii) ),
                                      Chen    = c(BFGS_estimates.Y.chen, unlist(model.Y.chen) ),
                                      Weibull = c(BFGS_estimates.Y.weibull, unlist(model.Y.weibull) )
                                      )
colnames(data.Y.results)[1:3]    <- c("alpha2", "r2", "s2")
round(data.Y.results, 5)
round(rbind(test.Y.burrxii, test.Y.chen, test.Y.weibull), 5)


##### Reliability estimations #####
# MLE 
N                                <- 50000 # sample size for the Monte Carlo integration

R_BFGS_estimate.B.C              <- R_MO(g1="burrxii", g2="chen", alpha1=BFGS_estimates.X.burrxii[1], r1=BFGS_estimates.X.burrxii[2], s1=BFGS_estimates.X.burrxii[3],
                                      alpha2=BFGS_estimates.Y.chen[1], r2=BFGS_estimates.Y.chen[2], s2=BFGS_estimates.Y.chen[3], N, previous.sample=FALSE, seed=1)$value
R_BFGS_estimate.B.C

R_BFGS_estimate.B.W              <- R_MO(g1="burrxii", g2="weibull", alpha1=BFGS_estimates.X.burrxii[1], r1=BFGS_estimates.X.burrxii[2], s1=BFGS_estimates.X.burrxii[3],
                                       alpha2=BFGS_estimates.Y.weibull[1], r2=BFGS_estimates.Y.weibull[2], s2=BFGS_estimates.Y.weibull[3], N, previous.sample=FALSE, seed=1)$value
R_BFGS_estimate.B.W


# Bayes
Non.informative.prior           <- rep(0.001,6)
Prior1.X                        <- rep(1,6)
Prior1.Y                        <- rep(1,6)
seed_number                     <- 12024
burnin.size                     <- 5000
mcmc.size                       <- 20000
thin.size                       <- 20

g1                              <- "burrxii" 
g2                              <- "chen"
Bayes.B.C.Prior0                <- Bayes.real.data(g1, g2, X, Y, BFGS_estimates.X = BFGS_estimates.X.burrxii, BFGS_estimates.Y = BFGS_estimates.Y.chen,
                                                   Prior.X = Non.informative.prior, Prior.Y = Non.informative.prior, burnin.size = 5000, mcmc.size = 20000, thin.size = 20, seed_number)

Bayes.B.C.PriorI                <- Bayes.real.data(g1, g2, X, Y, BFGS_estimates.X = BFGS_estimates.X.burrxii, BFGS_estimates.Y = BFGS_estimates.Y.chen,
                                                   Prior.X = Prior1.X, Prior.Y=Prior1.Y, burnin.size = 5000, mcmc.size = 20000, thin.size = 20, seed_number)

BC.estimates                    <- c( MLE = R_BFGS_estimate.B.C, B_Prior0 = Bayes.B.C.Prior0$R_estimate, B_PriorI = Bayes.B.C.PriorI$R_estimate)
round(BC.estimates, 5)


g1                              <- "burrxii" 
g2                              <- "weibull"
Bayes.B.W.Prior0                <- Bayes.real.data(g1, g2, X, Y, BFGS_estimates.X = BFGS_estimates.X.burrxii, BFGS_estimates.Y = BFGS_estimates.Y.weibull,
                                                   Prior.X = Non.informative.prior, Prior.Y = Non.informative.prior, burnin.size = 5000, mcmc.size = 20000, thin.size = 20, seed_number)

Bayes.B.W.PriorI                <- Bayes.real.data(g1, g2, X, Y, BFGS_estimates.X = BFGS_estimates.X.burrxii, BFGS_estimates.Y = BFGS_estimates.Y.weibull,
                                                   Prior.X = Prior1.X, Prior.Y = Prior1.Y, burnin.size = 5000, mcmc.size = 20000, thin.size = 20, seed_number)

BW.estimates                    <- c(MLE = R_BFGS_estimate.B.W, B_Prior0 = Bayes.B.W.Prior0$R_estimate, B_PriorI = Bayes.B.W.PriorI$R_estimate)
round(BW.estimates, 5)



