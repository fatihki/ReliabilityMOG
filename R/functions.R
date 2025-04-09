#================================================================================================================
# Functions in the application
#
# author: Fatih Kızılaslan (fkizilaslan@yahoo.com) and Gülce Cüran
# date: 22-December-2024
#================================================================================================================

# Reliability of MO family distributions calculation
R_MO <- function(g1, g2, alpha1, r1, s1, alpha2, r2, s2, N, previous.sample=FALSE, seed ){
  # previous.sample: TRUE/FALSE. You can specify the samples in Monte Carlo integration of R.
  # If you have a sample before, use TRUE otherwise FALSE.
  # seed: seed for the random number generation for Monte Carlo integration of R.
  
  if(g1=="weibull"){
    den=function(par,x){r=par[2]; s=par[3]; dweibull(x,shape=r,scale=s)}
    cum=function(par,x){r=par[2]; s=par[3]; pweibull(x,shape=r,scale=s)}
  }
  
  if(g1=="burrxii"){
    den=function(par,x){r=par[2]; s=par[3]; return(r*s*(x^(s-1))*(1+(x^s))^(-r-1))}
    cum=function(par,x){r=par[2]; s=par[3]; return(1-(1+(x^s))^(-r))}
  }
  
  if(g1=="chen"){
    den=function(par,x){r=par[2]; s=par[3]; return(r*s*(x^(r-1))*exp(x^r)*exp(s*(1-exp(x^r))))}
    cum=function(par,x){r=par[2]; s=par[3]; return(1-exp(s*(1-exp(x^r)))) }
  }
  
  f1_pdf<-function(par,x){
    f0=den(par,x)
    c0=cum(par,x)
    alpha=par[1]
    f = (alpha*f0)/((1-(1-alpha)*(1-c0))^2)
    return(f) }
  
  F1_cdf<-function(par,x){
    f0=den(par,x)
    c0=cum(par,x)
    alpha=par[1]
    FF = c0/(1-(1-alpha)*(1-c0))
    return(FF) }
  
  if(g2=="weibull"){
    den2=function(par,x){r=par[2]; s=par[3]; dweibull(x,shape=r,scale=s)}
    cum2=function(par,x){r=par[2]; s=par[3]; pweibull(x,shape=r,scale=s)}
  }
  
  if(g2=="burrxii"){
    den2=function(par,x){r=par[2]; s=par[3]; r*s*(x^(s-1))*(1+(x^s))^(-r-1)}
    cum2=function(par,x){r=par[2]; s=par[3]; 1-(1+(x^s))^(-r)}
  }
  
  if(g2=="chen"){
    den2=function(par,x){r=par[2]; s=par[3]; r*s*x^(r-1)*exp(x^r)*exp(s*(1-exp(x^r)))}
    cum2=function(par,x){r=par[2]; s=par[3]; 1-exp(s*(1-exp(x^r)))}
  }
  
  f2_pdf<-function(par,x){
    f0=den2(par,x)
    c0=cum2(par,x)
    alpha=par[1]
    g=(alpha*f0)/((1-(1-alpha)*(1-c0))^2)
    return(g) }
  
  F2_cdf<-function(par,x){
    f0=den2(par,x)
    c0=cum2(par,x)
    alpha=par[1]
    FF = c0/(1-(1-alpha)*(1-c0))
    return(FF) }
  
  r_integrand <-function(x){
    return( (1-F1_cdf(par=c(alpha1,r1,s1),x)) * f2_pdf(par = c(alpha2,r2,s2), x) ) }
  
  I = try( integrate(r_integrand, lower=0, upper=Inf), silent = TRUE)
  # Monte Carlo integration of R 
  
  # random sample from BurType XII dist.
  rburrxii <- function(n,r,s){
    u = runif(n)
    return( ( (1-u)^(-1/r) - 1 )^(1/s) ) }
  # random sample from Chen dist.
  rchen   <- function(n,r,s){
    u = runif(n)
    return( ( log( 1-(log(1-u))/s ) ) ^(1/r) ) }
  # H is a remaining part of the integrand of R except g2(y) pdf function
    H           <- function(y,alpha1,r1,s1,alpha2,r2,s2) {
    numerator   <- (1-F1_cdf(par=c(alpha1,r1,s1),x=y))*alpha2
    denominator <- ( 1-(1-alpha2)*(1-cum2(par = c(alpha2,r2,s2), x=y)) )^2
    return( numerator/denominator) }
  
  if ( class( I ) == "try-error") {
    if (  previous.sample == FALSE ){
      if(g2 == "weibull") {set.seed(seed); Y  = rweibull(N,shape=r2,scale=s2) }
      if(g2 == "burrxii") {set.seed(seed); Y  = rburrxii(N,r2,s2) }
      if(g2 == "chen")    {set.seed(seed); Y  = rchen(N,r2,s2) }
    } else {
      Y = previous.sample
    }

    MC_I           <- H(y=Y,alpha1,r1,s1,alpha2,r2,s2)           # values of H(y) wrt Y sample from g2 pdf
    S2             <- sum( (MC_I- mean(MC_I))^2) / (length(Y)-1) # sample variance of MC integral since it is approximation
    
    MC_integral    <- mean(MC_I)
    MC_integral_sample_variance <- S2
    
    out <- list(value = MC_integral, sample_var = MC_integral_sample_variance, sample=Y)
  } else{
    out <- list(value = I$value, sample = NA)
  }
  return(out)
}

MO_pdf      <- function(par, x, g){
  if (par[1]<0) return(-Inf)
  if (par[2]<0) return(-Inf)
  if (par[3]<0) return(-Inf)
  
  if(g == "weibull"){
    den = function(par,x){r=par[2]; s=par[3]; dweibull(x,shape=r,scale=s)}
    cum = function(par,x){r=par[2]; s=par[3]; pweibull(x,shape=r,scale=s)}
  }
  
  if(g == "burrxii"){
    den = function(par,x){r=par[2]; s=par[3]; return(r*s*(x^(s-1))*(1+(x^s))^(-r-1))}
    cum = function(par,x){r=par[2]; s=par[3]; return(1-(1+(x^s))^(-r))}
  }
  
  if(g=="chen"){
    den=function(par,x){r=par[2]; s=par[3]; return(r*s*(x^(r-1))*exp(x^r)*exp(s*(1-exp(x^r))))}
    cum=function(par,x){r=par[2]; s=par[3]; return(1-exp(s*(1-exp(x^r)))) }
  }
  
  f0       <- den(par,x)
  c0       <- cum(par,x)
  alpha    <- par[1]
  gpdf     <- (alpha*f0)/((1-(1-alpha)*(1-c0))^2)
  return(gpdf)
}


MO_cdf      <- function(x, par, g){
  if (par[1]<0) return(-Inf)
  if (par[2]<0) return(-Inf)
  if (par[3]<0) return(-Inf)
  
  if(g == "weibull"){
    den = function(par,x){r=par[2]; s=par[3]; dweibull(x,shape=r,scale=s)}
    cum = function(par,x){r=par[2]; s=par[3]; pweibull(x,shape=r,scale=s)}
  }
  
  if(g == "burrxii"){
    den = function(par,x){r=par[2]; s=par[3]; return(r*s*(x^(s-1))*(1+(x^s))^(-r-1))}
    cum = function(par,x){r=par[2]; s=par[3]; return(1-(1+(x^s))^(-r))}
  }
  
  if(g=="chen"){
    den=function(par,x){r=par[2]; s=par[3]; return(r*s*(x^(r-1))*exp(x^r)*exp(s*(1-exp(x^r))))}
    cum=function(par,x){r=par[2]; s=par[3]; return(1-exp(s*(1-exp(x^r)))) }
  }
  
  c0       <- cum(par,x)
  alpha    <- par[1]
  gcdf     <- c0 / (1-(1-alpha)*(1-c0))
  return(gcdf)
}


# negative log-likelihood
log_likelihood <- function(par, data, g){
  if (par[1]<0) return(-Inf)
  if (par[2]<0) return(-Inf)
  if (par[3]<0) return(-Inf)
  L <- -sum(log( MO_pdf(par, x=data, g)  ))
  if (is.na(L)==TRUE) {return(-Inf)} else {return(L)}
  return(L)
}
# gradient of negative log-likelihood
grad_log_likelihood <- function(par, data, g){
  return ( grad(func=log_likelihood, x=par, method="simple", data=data, g=g) )
}

# log-posterior
log_posterior <- function(par, data, hyperparameters) {
  g        <- hyperparameters$g
  alpha    <- par[1] # alpha
  r        <- par[2] # shape parameter
  s        <- par[3] # scale parameter
  if (par[1]<0) return(-Inf)
  if (par[2]<0) return(-Inf)
  if (par[3]<0) return(-Inf)
  a1<- hyperparameters$prior[1]; b1<- hyperparameters$prior[2];
  a2<- hyperparameters$prior[3]; b2<- hyperparameters$prior[4];
  a3<- hyperparameters$prior[5]; b3<- hyperparameters$prior[6];
  log.lik   <- sum(log(MO_pdf(par, x=data, g)))
  log.prior <- dgamma(alpha,a1,b1,log = TRUE) + dgamma(r,a2,b2,log = TRUE) + dgamma(s,a3,b3,log = TRUE) 
  L         <- log.lik + log.prior
  if (is.na(L) ==TRUE ) { return(-Inf)
  } else {
    return(L) }
}

# RMO_bootstrap function for the real data case
RMO_bootstrap_R1           <- function(g1, g2, n, m, alpha1.hat, r1.hat, s1.hat, alpha2.hat, r2.hat, s2.hat, R_BFGS.hat, B, B1, parallel = TRUE, seed ){
  
  R_bootstrap_BFGS         <- sd_Rstar_boot <- t_star <- c()
  
  if(parallel){
    boot.results           <- foreach( i=1:B, .combine=rbind) %dopar% {
      boot.fit             <- RMO_bootstrap_step_R1(g1, g2,  n, m, alpha1.hat, r1.hat, s1.hat, alpha2.hat, r2.hat, s2.hat, seed + i )
      # boot.t calculation
      boot.t.fit           <- lapply( 1:B1, function(j) { 
        RMO_bootstrap_step(g1, g2,  n, m, alpha1.hat = boot.fit$optim.X_BFGS_result$par[1], r1.hat = boot.fit$optim.X_BFGS_result$par[2], s1.hat = boot.fit$optim.X_BFGS_result$par[3],
                              alpha2.hat = boot.fit$optim.Y_BFGS_result$par[1], r2.hat = boot.fit$optim.Y_BFGS_result$par[2], s2.hat = boot.fit$optim.Y_BFGS_result$par[3], seed + j  )$R_bootstrap_BFGS
      })
      return( c( R_bootstrap = boot.fit$R_bootstrap_BFGS, sd_Rstar_boot = sd(unlist(boot.t.fit)) ) )
    }
    
    R_bootstrap_BFGS       <- as.numeric(boot.results[, 1])
    sd_Rstar_boot          <- as.numeric(boot.results[, 2])
  }
  
  t_star                   <- ( R_bootstrap_BFGS -  R_BFGS.hat ) /sd_Rstar_boot
  
  return( list( R_bootstrap_BFGS = R_bootstrap_BFGS, t_star = t_star ) )
}

# RMO_bootstrap_step_R1 for the real data case
RMO_bootstrap_step_R1     <- function(g1, g2,  n, m, alpha1.hat, r1.hat, s1.hat, alpha2.hat, r2.hat, s2.hat, seed ){
  
  set.seed(seed) 
  x.star                  <- rmog(n, spec=g1, beta=alpha1.hat, r=r1.hat, s = s1.hat)
  set.seed(seed + 2024) 
  y.star                  <- rmog(m, spec=g2, beta=alpha2.hat, r=r2.hat, s = s2.hat)
  
  initials.x.bootstrap    <- runif(1,0.75,1.25)*c(alpha1.hat, r1.hat, s1.hat)
  initials.y.bootstrap    <- runif(1,0.75,1.25)*c(alpha2.hat, r2.hat, s2.hat)
  
  # using optim with "BFGS" algorithm
  optim.X_BFGS_result     <- try( optim(par = initials.x.bootstrap, fn = log_likelihood, gr = grad_log_likelihood,
                                      g = g1, data = x.star, method = "BFGS" , control = list(maxit = 20000), hessian = F) )
  
  if ( class( optim.X_BFGS_result ) == "try-error") {
    set.seed(seed + 5 ) 
    x.star                <- rmog(n, spec=g1, beta=alpha1.hat, r=r1.hat, s = s1.hat)
    initials.x.bootstrap  <- runif(1,0.75,1.25)*c(alpha1.hat, r1.hat, s1.hat)
    optim.X_BFGS_result   <- try( optim(par = initials.x.bootstrap, fn = log_likelihood, gr = grad_log_likelihood,
                                        g = g1, data = x.star, method = "BFGS" , control = list(maxit = 20000), hessian = F) )
  }
  
  if ( class( optim.X_BFGS_result ) == "try-error") {
    set.seed(seed + 15 ) 
    x.star                <- rmog(n, spec=g1, beta=alpha1.hat, r=r1.hat, s = s1.hat)
    initials.x.bootstrap  <- runif(1,0.75,1.25)*c(alpha1.hat, r1.hat, s1.hat)
    optim.X_BFGS_result   <- try( optim(par = initials.x.bootstrap, fn = log_likelihood, gr = grad_log_likelihood,
                                        g = g1, data = x.star, method = "BFGS" , control = list(maxit = 20000), hessian = F) )
  }
  
  if ( class( optim.X_BFGS_result ) == "try-error") {
    next
  }
  
  optim.Y_BFGS_result    <- try( optim(par = initials.y.bootstrap, fn = log_likelihood, gr = grad_log_likelihood,
                                      g = g2, data = y.star, method = "BFGS" , control = list(maxit = 20000), hessian = F) )
  
  if ( class( optim.Y_BFGS_result ) == "try-error") {
    set.seed(seed + 2024 + 5) 
    y.star                <- rmog(m, spec=g2, beta=alpha2.hat, r=r2.hat, s = s2.hat)
    initials.y.bootstrap  <- runif(1,0.75,1.25)*c(alpha2.hat, r2.hat, s2.hat)
    optim.Y_BFGS_result   <- try( optim(par = initials.y.bootstrap, fn = log_likelihood, gr = grad_log_likelihood,
                                        g = g2, data = y.star, method = "BFGS" , control = list(maxit = 20000), hessian = F) )
  }
  
  if ( class( optim.Y_BFGS_result ) == "try-error") {
    set.seed(seed + 2024 + 15) 
    y.star                <- rmog(m, spec=g2, beta=alpha2.hat, r=r2.hat, s = s2.hat)
    initials.y.bootstrap  <- runif(1,0.75,1.25)*c(alpha2.hat, r2.hat, s2.hat)
    optim.Y_BFGS_result   <- try( optim(par = initials.y.bootstrap, fn = log_likelihood, gr = grad_log_likelihood,
                                        g = g2, data = y.star, method = "BFGS" , control = list(maxit = 20000), hessian = F) )
  }
  
  if ( class( optim.Y_BFGS_result ) == "try-error") {
    next
  }
  
  R_bootstrap_BFGS     <- R_MO(g1, g2, alpha1=optim.X_BFGS_result$par[1], r1=optim.X_BFGS_result$par[2], s1=optim.X_BFGS_result$par[3],
                               alpha2=optim.Y_BFGS_result$par[1], r2=optim.Y_BFGS_result$par[2], s2=optim.Y_BFGS_result$par[3], N, previous.sample=NA, seed=1)$value  #usage seed = 1 for the real data, in other case usage "seed"
  
  return( list( R_bootstrap_BFGS = R_bootstrap_BFGS,  optim.X_BFGS_result=optim.X_BFGS_result,  optim.Y_BFGS_result=optim.Y_BFGS_result ) )
  
}

# for real data analysis
Model_evaluation_MOG <- function(g, par, data) {
  # R^2:coefficient of determination, Root mean square (RMSE), K-S test values using i/(n+1), Akaike information criteria (AIC), Bayesian information criteria (BIC)
  # R^2, RMSE, KS, -2LogLik, AIC, BIC
  # par = c( alpha, r, s)
  sorted_data <- sort(data)
  w1          <- w2        <- w3       <- c()
  mean_estF   <- mean( MO_cdf(sorted_data, par, g ) )
  
  for (i in 1:length(data) ) {
    w1[i]     <- abs( MO_cdf(sorted_data[i], par, g ) - (i/(length(data)+1 )) )
    w2[i]     <- w1[i]^2                 # for R^2 value
    w3[i]     <- (w1[i] - mean_estF)^2   # for R^2 value
  }
  
  R2     <- 1- (sum(w2)/sum(w3))
  RMSE   <- sqrt(mean(w2))
  KS     <- max(w1)
  LL     <- log_likelihood(par, data, g) # negative log-likelihood
  AIC    <- (2*3)+(2*LL)
  BIC    <- (2*log(length(data))) + (2*LL)
  
  return( list( R2 = R2, RMSE = RMSE, KS = KS, NegLL = LL, AIC = AIC, BIC = BIC ) )
  
}

# function to run 3 MCMC chains given Priors
Bayes.real.data          <- function(g1, g2, X, Y, BFGS_estimates.X, BFGS_estimates.Y, Prior.X, Prior.Y, burnin.size=5000, mcmc.size=20000, thin.size=20, seed_number){
  
  posterior.samples1.X   <- try( MCMCmetrop1R( fun=log_posterior, theta.init = c(BFGS_estimates.X[1], BFGS_estimates.X[2], BFGS_estimates.X[3] ) ,
                                              data=X, hyperparameters=list(prior=Prior.X, g=g1),
                                              burnin=burnin.size, mcmc=mcmc.size, thin=thin.size, logfun=T,
                                              verbose=0, tune = 1, seed= seed_number ) )
  
  posterior.samples2.X  <- try( MCMCmetrop1R( fun=log_posterior, theta.init=c(BFGS_estimates.X[1], BFGS_estimates.X[2], BFGS_estimates.X[3] ) ,
                                              data=X, hyperparameters=list(prior=Prior.X, g=g1),
                                              burnin=burnin.size, mcmc=mcmc.size, thin=thin.size, logfun=T,
                                              verbose=0, tune = 1, seed= seed_number + 100 ) )
  
  posterior.samples3.X  <- try( MCMCmetrop1R( fun=log_posterior, theta.init=c(BFGS_estimates.X[1], BFGS_estimates.X[2], BFGS_estimates.X[3] ) ,
                                              data=X, hyperparameters=list(prior=Prior.X, g=g1),
                                              burnin=burnin.size, mcmc=mcmc.size, thin=thin.size, logfun=T,
                                              verbose=0, tune = 1, seed= seed_number + 200 ) )
  
  MH_estimates.X        <- rbind( as.matrix(posterior.samples1.X), as.matrix(posterior.samples2.X), as.matrix(posterior.samples3.X) )
  
  
  posterior.samples1.Y  <- try( MCMCmetrop1R(fun=log_posterior, theta.init= c(BFGS_estimates.Y[1], BFGS_estimates.Y[2], BFGS_estimates.Y[3]),
                                             data=Y,  hyperparameters=list(prior=Prior.Y, g=g2),
                                             burnin=burnin.size, mcmc=mcmc.size, thin=thin.size, logfun=T,
                                             verbose=0, tune = 1, seed= seed_number ) )
  
  posterior.samples2.Y  <- try( MCMCmetrop1R(fun=log_posterior, theta.init= c(BFGS_estimates.Y[1], BFGS_estimates.Y[2], BFGS_estimates.Y[3]),
                                             data=Y,  hyperparameters=list(prior=Prior.Y, g=g2),
                                             burnin=burnin.size, mcmc=mcmc.size, thin=thin.size, logfun=T,
                                             verbose=0, tune = 1, seed= seed_number + 100 ) )
  
  posterior.samples3.Y  <- try( MCMCmetrop1R(fun=log_posterior, theta.init= c(BFGS_estimates.Y[1], BFGS_estimates.Y[2], BFGS_estimates.Y[3]),
                                             data=Y,  hyperparameters=list(prior=Prior.Y, g=g2),
                                             burnin=burnin.size, mcmc=mcmc.size, thin=thin.size, logfun=T,
                                             verbose=0, tune = 1, seed= seed_number + 200 ) )
  
  MH_estimates.Y         <- rbind( as.matrix(posterior.samples1.Y), as.matrix(posterior.samples2.Y), as.matrix(posterior.samples3.Y) )
  
  
  Ri_MH_estimates          <- c()
  for (ij in 1:nrow(MH_estimates.X) ) {
    # R values based on the MCMC samples of the parameters
    Ri_MH_estimates[ij]   <- R_MO(g1, g2,alpha1=MH_estimates.X[ij,1], r1=MH_estimates.X[ij,2], s1=MH_estimates.X[ij,3],
                                alpha2=MH_estimates.Y[ij,1], r2=MH_estimates.Y[ij,2], s2=MH_estimates.Y[ij,3], N, previous.sample=FALSE, seed=1 )$value
  }
  
  HPD_R                   <- HPDinterval(as.mcmc(Ri_MH_estimates), prob = 0.95)
  R_HPD_lower             <- HPD_R[1]
  R_HPD_upper             <- HPD_R[2]
  R_estimate              <- mean(Ri_MH_estimates)  
  
  R_MH_chains             <- mcmc.list( mcmc( matrix( Ri_MH_estimates[1:(length(Ri_MH_estimates)/3)], ncol = 1, dimnames = list(NULL, "R")  ) ),
                                       mcmc( matrix( Ri_MH_estimates[((length(Ri_MH_estimates)/3)+1):(length(Ri_MH_estimates)*(2/3))] , ncol = 1, dimnames = list(NULL, "R") )),
                                       mcmc( matrix( Ri_MH_estimates[(length(Ri_MH_estimates)*(2/3)+1):length(Ri_MH_estimates) ] , ncol = 1, dimnames = list(NULL, "R") ))
                                       )
  R_GelmanF_estimate     <- gelman.diag(R_MH_chains , confidence = 0.95, autoburnin = F)$psrf[1] 
  
  return( list( R_estimate = R_estimate,  R_HPD = c(R_HPD_lower, R_HPD_upper),
                R_MH_chains.list = R_MH_chains, Rhat = R_GelmanF_estimate ) )
  
}
