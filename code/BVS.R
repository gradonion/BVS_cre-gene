#' @param y length N vec
#' @param X NxP mat
#' @param A lenght P vec
sample_gb = function(y,X,A,alpha,sigma2_e,sigma2_b,beta){
  
  N = dim(X)[1]
  P = dim(X)[2]
  Pgamma1 = rep(0,P)
  gamma = rep(0,P)
  # beta = rep(0,P)
  for (j in 1:P){
    lambda_j = 1/(sum(X[,j]^2)/sigma2_e + 1/sigma2_b)
    sum_xy = 0
    for (i in 1:N){
      sum_xy = sum_xy + X[i,j]*(y[i]-sum(X[i,]*beta)+X[i,j]*beta[j])
    }
    nu_j = lambda_j * sum_xy/sigma2_e
    Pgamma1[j] = sqrt(lambda_j/sigma2_b)*exp(nu_j^2/lambda_j/2 + alpha[1]*A[j]+alpha[2])
    Pgamma1[j] = 1-1/(Pgamma1[j]+1)
    gamma[j] = rbinom(1,1,Pgamma1[j])
    beta[j] = gamma[j] * rnorm(1,mean = nu_j,sd = sqrt(lambda_j))
  }
  return(list(gamma =  gamma, beta = beta, pi = Pgamma1))
}

#' @param prior_sigma2b sigma2_b ~ IG(g1,h1)
#' @param gamma length P
#' @param beta length P
#' @return sigma2_b: scalar
sample_sigma2_b = function(prior_sigma2b,gamma,beta){
  u = prior_sigma2b$g1 + sum(gamma)/2
  v = prior_sigma2b$h1 + sum(beta^2)/2
  sigma2_b = rigamma(1,u,v)
  return(sigma2_b)
}

#' @param prior_sigma2e sigma2_e ~ IG(g0,h0)
#' @param y length N
#' @param beta length P
#' @param X NxP matrix
#' @return sigma2_e: scalar
sample_sigma2_e = function(prior_sigma2e,y,beta,X){
  N = dim(X)[1]
  P = dim(X)[2]
  s = prior_sigma2e$g0 + N/2
  t = prior_sigma2e$h0 + sum( (y - t(X%*%beta)[1,])^2 )/2
  sigma2_e = rigamma(1,s,t)
  return(sigma2_e)
}

## EM for updating alpha (equivalent to logistic regression)
llh.alpha = function(alpha,x,y){
  return(-1*(alpha[1]*sum(y*x) + alpha[2]*sum(y) - sum(log(1+exp(alpha[1]*x+alpha[2])))))
}

## Complete likelihood 
LLH = function(y,X,A,alpha,sigma2_e,sigma2_b,beta,gamma){
  N = dim(X)[1]
  P = dim(X)[2]
  llh = -1 * (N/2*log(sigma2_e) +
                sum((y-(beta%*%t(X))[1,])^2)/(2*sigma2_e) +
                sum(gamma)*log(2*pi*sigma2_b)/2 +
                sum(gamma*(beta^2/sigma2_b/2)) +
                sum(gamma*log(1+exp(-alpha[1]*A-alpha[2]))+(1-gamma)*log(1+exp(alpha[1]*A+alpha[2]))) 
  )
  return(llh)
}

## Bayesian variable selction (BVS) implemented using Gibbs sampling:
#' @param y,X,A simulated data
#' @param alpha_true true enrichment parameter; if NULL, alpha will be estimated
#' @param sigma2_e_true,sigma2_b_true,beta_true true parameters that generated the data
#' @param prior_sigma2e,prior_sigma2b prior parameters on sigma2_b, sigma2_e in the bvs model
#' @param niter number of total gibbs sampling steps
#' @param burn number of burn-ins
#' @param seed random seed set for the process
#' @return A \code{"BVS"} object with the following elements:
#' \item{beta}{A vector of length P giving the estimated effect sizes
#'   for all P covariates (averaged over all iterations after burn-in).}
#' \item{PIP}{A vector of length P giving the posterior inclusion probabilities 
#'   for all P covariates (averaged over all iterations after burn-in).}
#' \item{sigma2_b}{The estimated prior variance for beta
#'   (averaged over all iterations after burn-in).}
#' \item{sigma2_e}{The estimated variance for error
#'   (averaged over all iterations after burn-in).}
#' \item{LLH}{}
#' 
BVS = function(y,X,A,alpha_true=NULL,
               prior_sigma2e=list(g0=1,h0=1),prior_sigma2b=list(g1=1,h1=1),
               niter=1000,burn=200,seed=12345){
  
  require(pscl)
  require(mvtnorm)
  
  ## Initialization:
  set.seed(seed)
  N = dim(X)[1] 
  P = dim(X)[2]
  sigma2_e = 4
  sigma2_b = 4
  Sigma2_e = rep(NA,niter)
  Sigma2_b = rep(NA,niter)
  Sigma2_e[1] = sigma2_e
  Sigma2_b[1] = sigma2_b
  
  Beta_update = matrix(NA,nrow = niter, ncol = P)
  beta = rep(0,P)
  Beta_update[1,] = beta
  Pgamma1_update = matrix(NA,nrow = niter, ncol = P)
  Pgamma1 = rep(0,P)
  Pgamma1_update[1,] = Pgamma1
  LLH.total = rep(NA,niter)
  
  ## Gibbs sampling iterations:
  iter = 1
  if(is.null(alpha_true)){
    alpha_update = matrix(NA,nrow = niter, ncol = 2)
    alpha = c(0,0)
    alpha_update[1,] = alpha
    while (iter<niter) {
      old_beta = beta
      gb = sample_gb(y,X,A,alpha,sigma2_e,sigma2_b,beta)
      beta = gb$beta
      gamma = gb$gamma
      Pgamma1 = gb$pi
      sigma2_b = sample_sigma2_b(prior_sigma2b,gamma,beta)
      sigma2_e = sample_sigma2_e(prior_sigma2e,y,beta,X)
      alpha = optim(par = c(0,0),fn = llh.alpha, x=A, y=Pgamma1, method = "BFGS")$par
      iter = iter+1
      Beta_update[iter,] = beta
      Pgamma1_update[iter,] = Pgamma1
      Sigma2_e[iter] = sigma2_e
      Sigma2_b[iter] = sigma2_b
      alpha_update[iter,] = alpha
      LLH.total[iter] = LLH(y,X,A,alpha,sigma2_e,sigma2_b,beta,gamma)
    }
  } else {
    alpha = alpha_true
    while (iter<niter) {
      gb = sample_gb(y,X,A,alpha,sigma2_e,sigma2_b,beta)
      beta = gb$beta
      gamma = gb$gamma
      Pgamma1 = gb$pi
      sigma2_b = sample_sigma2_b(prior_sigma2b,gamma,beta)
      sigma2_e = sample_sigma2_e(prior_sigma2e,y,beta,X)
      # alpha = optim(par = c(0,0),fn = llh.alpha, x=A, y=Pgamma1, method = "BFGS")$par
      iter = iter+1
      Beta_update[iter,] = beta
      Pgamma1_update[iter,] = Pgamma1
      Sigma2_e[iter] = sigma2_e
      Sigma2_b[iter] = sigma2_b
      # alpha_update[iter,] = alpha
      LLH.total[iter] = LLH(y,X,A,alpha,sigma2_e,sigma2_b,beta,gamma)
    }
  }
  
  ## Estimate the means of parameters
  # LLH_true = LLH(y,X,A,alpha_true,sigma2_e_true,sigma2_b_true,beta_true,gamma_true)
  # LLH.mean = mean(LLH.total[(burn+1):niter])
  sigma2_b.mean = mean(Sigma2_b[(burn+1):niter])
  sigma2_e.mean = mean(Sigma2_e[(burn+1):niter])
  beta.mean = colMeans(Beta_update[(burn+1):niter,])
  PIP.avg = colSums(Pgamma1_update[(burn+1):niter,])/(niter-burn)
  if(is.null(alpha_true)){
    alpha.mean = colMeans(alpha_update[(burn+1):niter,])
    result = list(beta = beta.mean, PIP = PIP.avg, alpha = alpha.mean,
                  sigma2_b = sigma2_b.mean, sigma2_e = sigma2_e.mean, LLH = LLH.total)
  } else {
    result = list(beta = beta.mean, PIP = PIP.avg, sigma2_b = sigma2_b.mean, 
                  sigma2_e = sigma2_e.mean, LLH = LLH.total)
  }
    
  return(result)
}
