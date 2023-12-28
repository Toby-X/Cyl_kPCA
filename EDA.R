library(readxl)
library(tidyverse)
library(cdmTools)
library(Rfast)

loc = read_excel("cylinder/point.xlsx")
X.train = read_excel("cylinder/X.train.xlsx")
Vx.train = read_excel("cylinder/pointx.xlsx")
Vy.train = read_excel("cylinder/pointu2.xlsx")
P.train = read_excel("cylinder/pointp.xlsx")

mu = rowMeans(Vx.train)
U_c = Vx.train-mu
U_c = as.matrix(U_c)
Sig = U_c%*%t(U_c)
U = HeteroPCA(Sig,5,max.iter = 1000)

ggplot(loc)+
  geom_point(aes(V1,V2))

paK(as.matrix(Vx.train))
Vx.train.mat = as.matrix(Vx.train)
paK(Vx.train.mat)
psych::fa.parallel(Vx.train.mat,fm="wls")
Vx.train.svd = svd(Vx.train)
which(cumsum(Vx.train.svd$d)/sum(Vx.train.svd$d)>0.99)[1]
Vy.train.svd = svd(Vy.train)

library(elasticnet)
Vx.dist = dist(Vx.train.mat)
Vx.dist = as.matrix(Vx.dist)
k.gaussian = exp(-Vx.dist/100)

K.central <- function(K){
  n = nrow(K)
  N1vec = matrix(rep(1/n,n),ncol=1)
  m1 = matrix(rep(t(N1vec)%*%K,n),nrow = n,byrow = T)
  m2 = matrix(rep(K%*%N1vec,n),nrow = n,byrow = F)
  m3 = matrix(rep(m1%*%N1vec,n),nrow = n,byrow = F)
  K.new = K-m1-m2+m3
  return(K.new)
}

K = K.central(k.gaussian)
sigma.sum = psych::tr(K)
K.evd = eigen(K)
K.evd$proportion = K.evd$values/sigma.sum
K.evd$idx = 1:nrow(K.evd$vectors)
K.evd$explianed = cumsum(K.evd$proportion)

K.df = data.frame(K.evd$idx,K.evd$proportion,K.evd$explianed)
colnames(K.df) = c("idx","proportion","explained")
ggplot(K.df[1:200,],aes(idx,explained))+
  geom_point()+
  geom_line(lwd=.75)+
  geom_hline(yintercept = .95, col="blue", lty = 2)+
  ylab("Explained Proportion")+
  xlab("Number of Principal Components")

K.svd = RSpectra::eigs(K,10)
cumsum(K.svd$values/sigma.sum)
leading.vector = K.svd$vectors[,1:10]
which(cumsum(K.svd$d)/sigma.sum>0.99)[1]
#960 and that is too many, how about controlling on the same level of PCs?

res1 = mat.mult(N1,k.gaussian)
res2 = matrix(rep(Crossprod(N1vec, k.gaussian),nrow(k.gaussian)),nrow=nrow(k.gaussian),byrow=T)

res1 = mat.mult(k.gaussian,N1)
res2 = matrix(rep(mat.mult(k.gaussian, N1vec),nrow(k.gaussian)),nrow=nrow(k.gaussian),byrow=F)

tic()
N1%*%k.gaussian
toc()

nmf_eu <- function(X,K,eps=1e-6,maxit=1e3){
  # mu = min(X)-1
  # X = X-mu
  err = 10
  iter = 0
  n = nrow(X)
  p = ncol(X)
  U = matrix(rbeta(n*K,1,1),nrow=n)
  score = matrix(rbeta(K*p,1,1),nrow=K)
  while (err >= eps & iter <= maxit) {
    iter = iter + 1
    score_tmp = score*(t(U)%*%X)/(t(U)%*%U%*%score)
    U_tmp = U*(X%*%t(score))/(U%*%score%*%t(score))
    err1 = max(abs(U_tmp-U))
    err2 = max(abs(score-score_tmp))
    err = max(err1,err2)
    U = U_tmp
    score = score_tmp
  }
  # return(list(U=U,score=score,iter=iter,mu=mu))
  return(list(U=U,score=score,iter=iter))
}

Vx_nmf = nmf_eu(as.matrix(Vx.train),5)
X = matrix(rpois(500*100,5),nrow=500)
Vx_nmf = nmf_eu(X,5)
