library(pscl)
library(mvtnorm)
library(VennDiagram)

## Load BVS and Elastic net modeling functions:
source("BVS.R")
source("ElasticNet.R")

## Load the function that generates simulation data:
source("data_gen.R")

## Data Initialization

# load("ATAC-seq_PeakData/peak_counts.cqn.Rdata") 
# peak.mean = colSums(peak_counts.cqn)/dim(peak_counts.cqn)[1]
# peak.cov = cov(as.matrix(peak_counts.cqn),method = 'pearson')

N = 50
P = 50
sparse = 0.2
niter = 1000
burn = 200

sigma2_e_true = 1  
alpha_true = c(0.5,-2)
prior_sigma2e = list(g0=1,h0=1)
prior_sigma2b = list(g1=1,h1=1)

binary.A = F
X.sd = 1

Sigma2_b = c(0.2,seq(0.5,4,0.5))

stats.summary = data.frame(matrix(nrow = length(Sigma2_b), ncol = 9))
names(stats.summary) = c('sigma2_b','sparse','en.precision','en.recall','bvs.precision',
                         'bvs.recall','bvs.inv_spec','bvs.pp0','bvs.pp1')
rep = 30
seed0 = 12345
set.seed(seed0)
seeds = sample(100000,size=rep)

for (s in 1:length(Sigma2_b)){

  sigma2_b_true = Sigma2_b[s]
  print(paste("True sigma2_b:",sigma2_b_true))
  stats.summary$sigma2_b[s] = sigma2_b_true
  
  en = list(overlap = matrix(nrow = rep,ncol = 3))
  bvs = list(overlap = matrix(nrow = rep,ncol = 3), means = matrix(nrow = rep,ncol = 4), pip = matrix(nrow = rep,ncol = 2))

  for (i in 1:rep){
    seed = seeds[i]
    data = data.gen(N,P,peak.mean=NULL,peak.cov=NULL,X.sd,binary.A,sparse,alpha_true,sigma2_e_true,sigma2_b_true,seed)
    X = data$X
    y = data$y
    A = data$A
    beta_true = data$beta_true
    en$overlap[i,] = EN.analysis(X,y,beta_true,seed)
    stats = BVS(y,X,A,alpha_true,sigma2_e_true,sigma2_b_true,beta_true,prior_sigma2e,prior_sigma2b,niter,burn,seed)
    bvs$overlap[i,] = stats$overlap
    bvs$means[i,] = stats$means
    bvs$pip[i,] = stats$pip
  }
  
  stats.summary$sparse[s] = mean(bvs$overlap[,1])/P
  print(paste("Averaged sparsity:",mean(bvs$overlap[,1])/P))
  stats.summary$en.precision[s] = mean(en$overlap[,3]/en$overlap[,2])
  stats.summary$en.recall[s] = mean(en$overlap[,3]/en$overlap[,1])
  stats.summary$bvs.precision[s] = mean(bvs$overlap[,3]/bvs$overlap[,2])
  stats.summary$bvs.recall[s] = mean(bvs$overlap[,3]/bvs$overlap[,1])
  stats.summary$bvs.inv_spec[s] = mean((bvs$overlap[,2]-bvs$overlap[,3])/(P-bvs$overlap[,1]))
  stats.summary$bvs.pp0[s] = mean(bvs$pip[,1])
  stats.summary$bvs.pp1[s] = mean(bvs$pip[,2])
}

plot(stats.summary$sigma2_b, stats.summary$bvs.precision, ylim = c(0,1), type = 'l',
     xlab = 'sigma2_b', ylab = '', col = 'indianred1', main = paste('X.sd',X.sd), 
     cex.axis = 1.2)

points(stats.summary$sigma2_b, stats.summary$bvs.precision, pch=19, col = 'indianred1')
lines(stats.summary$sigma2_b, stats.summary$bvs.recall, type = 'l', col = 'indianred1')
points(stats.summary$sigma2_b, stats.summary$bvs.recall, pch = 17, col = 'indianred1')

lines(stats.summary$sigma2_b, stats.summary$en.precision, type = 'l',col='turquoise3')
points(stats.summary$sigma2_b, stats.summary$en.precision,pch=19, col = 'turquoise3')
lines(stats.summary$sigma2_b, stats.summary$en.recall, type = 'l',col='turquoise3')
points(stats.summary$sigma2_b, stats.summary$en.recall, pch = 17, col='turquoise3')

legend(x=2.5,y=0.3,legend = c('BVS precision','BVS recall','EN precision','EN recall'), 
       col = c('indianred1','indianred1','turquoise3','turquoise3'),
       pch = c(19,17,19,17), bty = "n")



