library(foreach)
library(doParallel)
library(betareg)
library(mc2d)
library(ggplot2)
library(stats)

#' Inverse logit in the flexible Beta regression model
#' 
#' @param    X matrix; first columns contains only 1's
#' @param beta numeric vector; length 2
#' 
mu_logit <- function(X, beta){
  exp( as.numeric(X%*%beta) )/(1 + exp(as.numeric(X%*%beta)))
}

#' Density of beta flexible random variable
#' 
#' @param           y vector
#' @param       omega vector same length of y
#' @param           p numeric in interval (0,1)
#' @param          mu_logit vector length(y)
#' @param         phi numeric positive
#' @param mixtureType numeric either 0 or 1; describes mixtures in the beta flexible model
#' 
dBetaFlexible <- function(y,omega,p,mu,phi,mixtureType=c(0,1)){
  out <- numeric(length(y))
  
  if(mixtureType==0){
    lambda <- mu - p * omega
  } else {
    lambda <- mu + (1-p) * omega
  }
  
  indice <- 1:length(y)
  indices1 <- indice[y==1]
  indices0 <- indice[y==0]
  indices <- c(indices0, indices1)
  
  if(length(indices)==0){
    out <- dbeta(x=y,shape1=lambda*phi,shape2=(1-lambda)*phi)
  } else {
    out[indices] <- 1e-4
    out[-c(indices)] <- dbeta(y[-c(indices)],shape1=lambda[-c(indices)]*phi,
                              shape2=(1-lambda[-c(indices)])*phi) 
  }
  
out  
}

#' Sampling from posterior distribution of auxiliar parameter \code{\eta}
#' 
#' @param     y vector
#' @param omega vector
#' @param     p numeric in (0,1)
#' @param    mu vector
#' @param   phi numeric positive
#' 
posterior_v <- function(y,omega,p,mu,phi){
  A <- p * dBetaFlexible(y=y,omega=omega,p=p,mu=mu,phi=phi,mixtureType=1)
  B <- (1-p) * dBetaFlexible(y=y,omega=omega,p=p,mu=mu,phi=phi,mixtureType=0)
  
  probVector <- A/(A+B)
  
  bernTrials <- rbern(n=length(y),prob=probVector)
  n0 <- sum(bernTrials==0)
  n1 <- sum(bernTrials==1)
  
list(n0 = n0, n1 = n1, trials = bernTrials)
}

#' Sampling from posterior distribution of parameter \code{p}
#' 
#' @param n0 numeric integer
#' @param n1 numeric integer
#' 
posterior_p <- function(n0,n1){
  rbeta(1, shape1 = n1+1, shape2 = n0+1)
}

#' Sampling from posterior distribution of parameter \code{\omega}
#' 
#' @param  p numeric in (0,1)
#' @param mu vector
#' 
posterior_omega <- function(p,mu){
  Omega <- matrix(nrow=length(mu), ncol = 2)
  Omega[,1] <- mu/p
  Omega[,2] <- (1-mu)/(1-p)
  minimos <- apply(Omega, 1, min)

runif(length(mu),min=1e-4,max=0.9999) * minimos
}

#' Performs thinning in sequence of pseudo random numbers aiming to remove
#' the existing serial dependence
#' 
#' @param        lag numeric integer; controls degree of existing dependence; how separated the samples should be?
#' @param      delta numeric integer; companion parameter for lag
#' @param sampleSize numeric integer; length of samples after thinning
#' @param  chainSize numeric integer; length of original sequence
#'   
thinning <- function(lag,delta,sampleSize,chainSize){
  size <- floor(chainSize/lag)
  
  out <- c()
  out[1] <- sample((lag-delta):(lag+delta),1)
  while(length(out) < sampleSize){
    for(i in 2:size){
      out <- c(out,sample((i*lag-delta):(i*lag+delta),1))
    }
  }
  
  sample <- sample(out,sampleSize)
  
  list(sample=sample,burns=out,size=size)
}

#' Sampling from posterior distribution of parameter \code{\phi}
#' 
#' @param            y vector
#' @param            v vector
#' @param            n numeric sample size
#' @param           mu vector
#' @param        omega vector
#' @param            p numeric in (0,1)
#' @param phi_proposal numeric positive
#' @param        stDev numeric positive; determines standard deviation of jumping 
#'                     distribution
#' @param            g numeric positive; rate parameter of prior gamma distribution
#' @param            k numeric integer; used as a factor for scale parameter of 
#'                     prior gamma distribution
#' 
getSamplePosteriorPhi <- function(y,v,mu,omega,p,phi_proposal,stDev,g,k){
  acceptance <- 0
  output <- numeric(1)
  y_s <- y[v$trials==1] #Rescatamos las que pertenecen a Ys (1)
  y_t <- y[v$trials==0] #Rescatamos las que pertenecen a Yt (0)
  mu_s <- mu[v$trials==1]
  mu_t <- mu[v$trials==0]
  omega_s <- omega[v$trials==1]
  omega_t <- omega[v$trials==0]
  
  phiNew <- phi_proposal + rnorm(1,0,sd=stDev)
  phis <- c(phi_proposal, phiNew)
  
  alpha_aux <- sapply(1:length(phis), 
                      function(s) getLogLikeCond(phi=phis[s],p=p,
                                                 y_s=y_s,y_t=y_t,
                                                 omega_s=omega_s,omega_t=omega_t,
                                                 mu_s=mu_s,mu_t=mu_t))
  
  alpha <- min(0, alpha_aux[2]-alpha_aux[1]+
                 log(dgamma((phis[2]),k*g,g))-log(dgamma((phis[1]),k*g,g))) # +
  
  u <- log(runif(1,0,1))
  if(u <= alpha){
    output <- phiNew
    acceptance <- 1
  }else{
    output <- phi_proposal
  }
list(sample=output,acceptance=acceptance)  
}

#' Conditional log-likelihood function of Beta Flexible regression model
#' 
#' @param     phi numeric positive
#' @param       p numeric in (0,1)
#' @param     y_s vector
#' @param     y_t vector
#' @param omega_s vector; length(y_s)
#' @param omega_t vector; length(y_t)
#' @param    mu_s vector; length(y_s)
#' @param    mu_t vector; length(y_t)
#' 
getLogLikeCond <- function(phi,p,y_s,y_t,omega_s,omega_t,mu_s,mu_t){
  
  f1 <- logBetaFlexible(y=y_s,omega=omega_s,p=p,mu=mu_s,phi=phi,mixtureType=1)
  f0 <- logBetaFlexible(y=y_t,omega=omega_t,p=p,mu=mu_t,phi=phi,mixtureType=0)
  
  if(length(y_s)==0){
    f1 <- 0
  } 
  
  if(length(y_t)==0){
    f0 <- 0
  }
  
  f1 + f0
}

#' Log-likelihood function of beta flexible random variable.
#' 
#' @param           y vector
#' @param       omega vector
#' @param           p numeric in (0,1)
#' @param          mu vector
#' @param         phi numeric positive
#' @param mixtureType numeric either 0 or 1; describes mixtures in the beta flexible model
#' 
logBetaFlexible <- function(y,omega,p,mu,phi,mixtureType){
  
  out <- numeric(length(y))
  
  if(mixtureType==1){
    lambda <- mu + (1-p)*omega
  } else {
    lambda <- mu - p*omega
  }
  
  indice_pos <- (1:length(y))[y==1]
  indice_zero <- (1:length(y))[y==0]
  indices <- c(indice_pos, indice_zero)
  
  if(length(indices)==0){
    out <- getLogLikeBeta(y=y,a=lambda,b=phi) # need to define log Beta
  } else {
    # rep(log(1e-4),times=length(indices)) +
    out <- log(1e-4)*length(indices) +
      getLogLikeBeta(y[-c(indices)],a=lambda[-c(indices)],b=phi)
  }
  
  out  
}

#' Log-likelihood function of beta flexible random variable with media-precision parametrization
#' 
#' @param y vector
#' @param a vector; length(y)
#' @param b vector; length(y)
#'  
getLogLikeBeta <- function(y,a,b){
  length(y)*log(gamma(b))-
    sum(log(gamma(a*b)))-
    sum(log(gamma((1-a)*b)))+
    sum((a*b-1)*log(y)+((1-a)*b-1)*log(1-y))
}


#' Sampling posterior distribution of parameter \code{\beta}
#' 
#' @param             x vector
#' @param             y vector 
#' @param             v vector; length(y)
#' @param             p numeric in (0,1)
#' @param         omega vector; length(y)
#' @param           phi numeric positive
#' @param            mu vector length(y)
#' @param beta_proposal vector 2x1; initial value of MCMC
#' @param     stDevJump numeric vector; length 2; standard deviation of jumping
#'                      distribution in M-H algorithm
#' @param    stDevPrior numeric vector; length 2; standard deviation of prior
#'                      distribution of \code{beta}
#' 
getSamplePosteriorBeta <- function(x,y,v,p,omega,phi,mu,beta_proposal,stDevJump,
                                   stDevPrior){
  acceptance <- 0
  output <- c(NA,NA)
  y_s <- y[v$trials==1]
  y_t <- y[v$trials==0] 
  x_s <- x[v$trials==1]
  x_t <- x[v$trials==0]
  X_s <- cbind(rep(1,length(x_s)),x_s)
  X_t <- cbind(rep(1,length(x_t)),x_t)
  mu_s <- mu[v$trials==1]
  mu_t <- mu[v$trials==0]
  omega_s <- omega[v$trials==1]
  omega_t <- omega[v$trials==0]
  X <- cbind(rep(1,length(x)),x)
  
  betaNew <- beta_proposal + t(rmvnorm(1,mean=rep(0,2),sigma=diag(stDevJump)))
  
  alphaOld <- getLogLikeCond(phi=phi,p=p,y_s=y_s,y_t=y_t,
                              omega_s=omega_s,omega_t=omega_t,
                              mu_s=mu_s,mu_t=mu_t) + 
      log(dmvnorm(t(beta_proposal),mean=rep(0,2),sigma=diag(stDevPrior)))
  
  muNew_s <- mu_logit(X=X_s,beta=betaNew)
  muNew_t <- mu_logit(X=X_t,beta=betaNew)
  omegaNew_s <- posterior_omega(p=p,mu=muNew_s)
  omegaNew_t <- posterior_omega(p=p,mu=muNew_t)
  
  alphaNew <- getLogLikeCond(phi=phi,p=p,y_s=y_s,y_t=y_t,
                              omega_s=omegaNew_s,omega_t=omegaNew_t,
                              mu_s=muNew_s,mu_t=muNew_t) + 
      log(dmvnorm(t(betaNew),mean=rep(0,2),sigma=diag(stDevPrior)))
  
  alpha <- min(0, alphaNew-alphaOld)
  u <- log(runif(1,0,1))
  if(u <= alpha){
    output <- betaNew
    acceptance <- acceptance + 1
  } else {
    output <- beta_proposal
  }
  
list(sample=output,acceptance=acceptance)  
}

# ---

#' MCMC sampling to estimate parameters in the Flexible Beta Regression model.
#' Sampling strategy is based on Gibbs sampling to iteratively get posterior
#' distributions of parameters \code{p}, \code{omega}, \code{phi} and \code{\beta}.
#' A Metropolis-Hastings algorithm is required to sample from (conditional) posterior distribution
#' of \code{phi} and \code{beta}. Yielded samples are based on runs of several Markov chains.
#' Parallel computing is used to get several MCMC at once.
#' 
#' @param              x numeric vector
#' @param              y numeric vector
#' @param          p_ini numeric in (0,1); starting point for sampling distribution 
#'                       of \code{p}
#' @param      omega_ini numeric vector; starting point for sampling distribution 
#'                       of \code{omega}
#' @param        phi_ini numeric non-negative; starting point for sampling distribution 
#'                       of \code{phi}
#' @param       beta_ini numeric vector; length 2; starting point for sampling 
#'                       distribution of \code{beta}
#' @param      numChains numeric integer; number of MCMC required to sample from 
#'                       distribution of \code{phi} and \code{beta}
#' @param samplesByChain numeric integer; size of samples produced by each MCMC
#' @param       stDevPhi numeric non-negative; standard deviation of jumping distribution 
#'                       used in MH to sample from distribution of \code{phi}
#' @param              g numeric non-negative; auxiliar parameter in MH to sample 
#'                       from distribution of \code{phi}
#' @param              k numeric integer; auxiliar parameter in MH to sample from 
#'                       distribution of \code{phi}
#' @param  stDevBetaJump numeric vector; length 2; standard deviation of jumping 
#'                       distribution used in MH to sample from distribution of 
#'                       \code{beta}
#' @param stDevBetaPrior numeric vector; length 2; standard deviation of prior
#'                       distribution of \code{beta} used in MH algorithm                        
#'        
multipleGibbsFBR <- function(x=NULL,y,p_ini,omega_ini,phi_ini,beta_ini,
                             numChains=5,samplesByChain=100,numCores,
                             stDevPhi,g,k,
                             stDevBetaJump,stDevBetaPrior){
  
  cluster <- parallel::makeCluster(numCores, outfile = "")
  registerDoParallel(cluster) 
  
  clusterExport(cluster, c('mu_logit','dBetaFlexible','posterior_v','posterior_p',
                           'posterior_omega',
                           'getLogLikeCond','logBetaFlexible','getLogLikeBeta',
                           'gibbsFBR','getSamplePosteriorPhi','getSamplePosteriorBeta'))
  
  output <- foreach(1:numChains,.packages=c('mc2d', 'stats', 'betareg')) %dopar% {
    
    TEMP <- gibbsFBR(x=x,y=y,p_ini=p_ini,omega_ini=omega_ini,
                             phi_ini=phi_ini,beta_ini=beta_ini,
                             samples=samplesByChain,
                             stDevPhi,g,k,stDevBetaJump,stDevBetaPrior)
    
    return(TEMP)
  }

  stopCluster(cluster)
  
output
}

gibbsFBR <- function(x,y,p_ini,omega_ini,phi_ini,beta_ini,
                     samples,
                     stDevPhi,g,k,
                     stDevBetaJump,stDevBetaPrior){
  p <- numeric(samples)
  omega <- matrix(nrow = samples, ncol = length(x))
  phi <- numeric(samples)
  beta <- matrix(nrow = 2, ncol = samples)
  mu <- matrix(nrow = samples, ncol = length(x))
  
  p[1] <- p_ini
  phi[1] <- phi_ini
  beta[,1] <- beta_ini
  X <- cbind(rep(1,length(x)),x)
  mu[1,] <- mu_logit(X=X,beta=beta[,1])
  Omega <- matrix(nrow = length(y), ncol = 2)
  Omega[,1] <- mu[1,]/p[1]
  Omega[,2] <- (1-mu[1,])/(1-p[1])
  omega[1,] <- omega_ini * apply(Omega,1,min)
  
  acceptancePhi <- numeric(samples)
  acceptancePhi[1] <- 1
  acceptanceBeta <- numeric(samples)
  acceptanceBeta[1] <- 1
    
  for(i in 2:samples){
    v_post <- posterior_v(y=y,omega=omega[i-1,],p=p[i-1],mu=mu[i-1,],phi=phi[i-1])
    p[i] <- posterior_p(n0=v_post$n0,n1=v_post$n1)
    mu[i,] <- mu_logit(X=X,beta=beta[,i-1]) # CAMBIO AHORA
    omega[i,] <- posterior_omega(p=p[i],mu=mu[i,])
    
    TEMP_phi <- getSamplePosteriorPhi(y=y,v=v_post,mu=mu[i,],omega=omega[i,],
                                      p=p[i],phi_proposal=phi[i-1],stDev=stDevPhi,g=g,k=k)
    
    phi[i] <- TEMP_phi$sample
    acceptancePhi[i] <- TEMP_phi$acceptance

    TEMP_beta <- getSamplePosteriorBeta(x=x,y=y,v=v_post,p=p[i],omega=omega[i,],
                                        phi=phi[i],mu=mu[i,],
                                        beta_proposal=beta[,i-1],
                                        stDevJump=stDevBetaJump,
                                        stDevPrior=stDevBetaPrior)
    
    beta[,i] <- TEMP_beta$sample #beta_TEMP[,sample(1:ncol(beta_TEMP),1)]
    acceptanceBeta[i] <- TEMP_beta$acceptance
  }
  
  acceptanceRatePhi <- sum(acceptancePhi)/samples
  acceptanceRateBeta <- sum(acceptanceBeta)/samples
  
  list(p=p,omega=omega,phi=phi,acceptancePhi=acceptanceRatePhi,
       beta0=beta[1,],beta1=beta[2,],acceptanceBeta=acceptanceRateBeta)
}