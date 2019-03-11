## Function to generate simulation data:
data.gen = function(N,P,peak.mean=NULL,peak.cov=NULL,X.sd=1,binary.A=FALSE,
                    sparse=0.2,alpha_true,sigma2_e_true=1,sigma2_b_true=1,seed=12345){
  
  set.seed(seed)
  if (binary.A == T){
    A = rep(0,P)
    A[sample(1:P,floor(P*sparse))] = 1
  } else {
    A =  rnbinom(P, size=20, prob=0.25)
    A[sample(1:P,floor(P*(1-sparse)))] = 0
    A = log(A+1)
  }
  
  pi_1 = 1/(exp(-alpha_true[1]*A-alpha_true[2])+1)
  gamma_true = rbinom(P,1,pi_1)
  # print(sum(gamma_true==1))
  if (is.null(peak.cov)){
    X = matrix(rnorm(N*P,mean = 0,sd = X.sd), nrow = N, ncol = P) 
  } else {
    print('Using sample mean and covariance computed from real peak count data.')
    X = t(rmvnorm(P,peak.mean,peak.cov))
    if (N != length(peak.mean)){
      print(paste('Adjusting sample mean and covariance to',N,'dimension...'))
      fold = floor(N/length(peak.mean))
      tmp = X
      for (f in 2:fold){
        tmp = rbind(tmp,t(rmvnorm(P,peak.mean,peak.cov)))
      }
      X = tmp
    }
  }
  
  beta_true = gamma_true*rnorm(P,mean = 0,sd = sqrt(sigma2_b_true))
  e = rnorm(N,mean = 0,sd = sqrt(sigma2_e_true))
  y = beta_true%*%t(X)
  y = y[1,] + e
  return(list(y = y, X = X, A = A, beta_true = beta_true))
}
