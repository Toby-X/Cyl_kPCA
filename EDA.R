library(readxl)
library(tidyverse)
library(cdmTools)
library(Rfast)

loc = read_excel("cylinder/point.xlsx")
X.train = read_excel("cylinder/X.train.xlsx")
Vx.train = read_excel("cylinder/pointx.xlsx")
Vy.train = read_excel("cylinder/pointu2.xlsx")
P.train = read_excel("cylinder/pointp.xlsx")

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
k.gaussian = exp(-Vx.dist)

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
K.svd = svd(K)
which(cumsum(K.svd$d)/sigma.sum>0.99)[1]
#960 and that is too many, how about controlling on the same level of PCs?

res1 = mat.mult(N1,k.gaussian)
res2 = matrix(rep(Crossprod(N1vec, k.gaussian),nrow(k.gaussian)),nrow=nrow(k.gaussian),byrow=T)

res1 = mat.mult(k.gaussian,N1)
res2 = matrix(rep(mat.mult(k.gaussian, N1vec),nrow(k.gaussian)),nrow=nrow(k.gaussian),byrow=F)

tic()
N1%*%k.gaussian
toc()
