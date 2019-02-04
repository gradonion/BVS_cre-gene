EN.analysis = function(X,y,beta_true,seed){
  require(glmnet)
  require(VennDiagram)
  set.seed(seed)
  fit.elnet.cv = cv.glmnet(X, y, type.measure="mse", alpha=.5, family="gaussian")
  # plot(fit.elnet.cv)
  
  lambda.1se = fit.elnet.cv$lambda.1se
  lambda.min = fit.elnet.cv$lambda.min
  
  fit.elnet = glmnet(X, y, family="gaussian", alpha=.5)
  # plot(fit.elnet, xvar="lambda")
  # abline(v = c(log(lambda.1se),log(lambda.min)), lty = 2 )
  
  cv_error = fit.elnet.cv$cvm[fit.elnet.cv$lambda==lambda.1se]
  cv_sd = fit.elnet.cv$cvsd[fit.elnet.cv$lambda==lambda.1se]
  fit.1se = glmnet(X, y, lambda = lambda.1se, alpha=.5, family="gaussian")
  
  nonzero_beta = fit.1se$beta[fit.1se$beta[,1]!=0, 1]
  nonzero.indx = (1:P)[fit.1se$beta[,1]!=0]
  true.indx = (1:P)[beta_true!=0]
  overlap = calculate.overlap( x = list(true.indx,nonzero.indx))
  stats = as.numeric(sapply(overlap, length))
  
  return(stats)
}

