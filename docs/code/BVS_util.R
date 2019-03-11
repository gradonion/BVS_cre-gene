BVS_plot.PIP = function(pip,thres_pip=NULL,beta_true,thres_b=NULL,pos=NULL,...){

  if(is.null(pos)){
    pos = 1:length(pip)
  } else { 
    beta_true = beta_true[pos]
    pip = pip[pos]
  }
  
  plot(pos,pip,col="black", xlab='covariate index', ylab='PIP', pch=16, ...)
  
  if(is.null(thres_pip)){
    lines(c(0,length(pip)+1),rep(0.5,2),col='coral1',lty=2,lwd=1.2)
  } else {
    lines(c(0,length(pip)+1),rep(thres_pip,2),col='coral1',lty=2,lwd=1.2)
  }
  
  if(is.null(thres_b)){
    points(pos[beta_true!=0],pip[beta_true!=0],col='red',pch=16)
  } else {
    points(pos[abs(beta_true)>thres_b],pip[abs(beta_true)>thres_b],col='red',pch=16)
  }
}

BVS_selection.stats = function(pip,thres_pip=NULL,beta_true,thres_b=NULL){
  require(VennDiagram)
  P = length(pip)
  if(is.null(thres_b)){
    true.indx = (1:P)[beta_true!=0]
  } else {
    true.indx = (1:P)[abs(beta_true)>thres_b]
  }
  
  if(is.null(thres_pip)){
    est.indx = (1:P)[pip>0.5]
  } else {
    est.indx = (1:P)[pip>thres_pip]
  }
  
  pp = c(mean(pip[-true.indx]), mean(pip[true.indx]))
  names(pp) = c('pp0','pp1')
  
  intersect = match(est.indx,true.indx)
  overlap = c(length(true.indx),length(est.indx),length(na.omit(intersect)))
  names(overlap) = c('true signals','estimated signals','overlap')
    
  return(list(pp = pp, overlap = overlap))
}
