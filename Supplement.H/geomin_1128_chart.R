
#produce tabl 1.1 and h.2 in appendix
library(GPArotation)
source('cl1121_Lp_fun.R')
Qt=matrix(c(sqrt(3)/3,sqrt(3)/3,-sqrt(3)/3,0,1,0,0,-sqrt(2)/2,sqrt(2)/2),3,3,byrow=T)
solve(Qt)
X = matrix(0, 15, 3)
set.seed(1)
diag(X[1:3,])=runif(3,0.5,1)
diag(X[4:6,])=runif(3,0.5,1)
diag(X[7:9,])=runif(3,0.5,1)
diag(X[10:12,])=runif(3,0.5,1)
diag(X[13:15,])=runif(3,0.5,1)

Dprime = matrix(runif(18,0.5,1), 6,3)
diag(Dprime[1:3,])=0
diag(Dprime[4:6,])=0

Aprime = rbind(X%*%solve(Qt),Dprime) # A=Aprime%*%Qt
A  = rbind(X, Dprime %*% Qt)
sum(Aprime==0)
sum(A==0)

geominf <- function(L,delta=0.01){
  sum(exp(rowSums(log(L^2 + delta))/ncol(L)))}
delta=10^(-6) #0.000001
geominf(Aprime,delta=delta)
geominf(A,delta=delta)

f.quartimin <- function(L){
  X <-  L^2 %*% (!diag(TRUE,ncol(L))) 
  return(sum(L^2 * X)/4)
}

f.simplimax <- function(L, k=nrow(L)){
  # k: Number of close to zero loadings
  Imat <- sign(L^2 <= sort(L^2)[k])
  return(sum(Imat*L^2)) 
}

A_geo=geominQ(A,maxit=100000)$loadings
A_L1 = irls(1,A)$L
res=list()

res$l1$mse= mean((irls(1,A)$L-A)^2)
res$l0.5$mse= mean((irls(0.5,A)$L-A)^2)
res$oblimin$mse= mean((oblimin(A)[[1]]-A)^2)
res$quartimin$mse= mean((quartimin(A)[[1]]-A)^2)
res$simplimax$mse= mean((simplimax(A)[[1]]-A)^2)
res$geominQ$mse= mean((A_geo-A)^2)
res$geomin0$mse= mean((Aprime-A)^2)
res$promax$mse=mean((promax(A)[[1]]-A)^2)

res$l1$L0= sum(abs(irls(1,A)$L)<0.01)
res$l0.5$L0= sum(abs(irls(0.5,A)$L)<0.01)
res$oblimin$L0= sum(abs(oblimin(A)[[1]])<0.01)
res$quartimin$L0=sum(abs(quartimin(A)[[1]])<0.01)
res$simplimax$L0= sum(abs(simplimax(A)[[1]])<0.01)
res$geominQ$L0= sum(abs(A_geo)<0.01)
res$geomin0$L0= sum(abs(Aprime)<0.01)
res$promax$L0=sum(abs(promax(A)[[1]])<0.01)

# the obj value at A, where the function was copied form GPArotation package for comparison
res$l1$fA= sum(abs(A))
res$l0.5$fA=sum(abs(A)^(0.5))
res$oblimin$fA= f.quartimin(A)
res$quartimin$fA=f.quartimin(A)
res$simplimax$fA= f.simplimax(A)
res$geominQ$fA= geominf(A)
res$geomin0$fA= geominf(A,delta=0)
res$promax$fA=NA

# the obj value at rot(A)
res$l1$f= sum(abs(irls(1,A)$L))
res$l0.5$f=sum(abs(irls(0.5,A)$L)^(0.5))
res$oblimin$f= f.quartimin(oblimin(A)[[1]])
res$quartimin$f=f.quartimin(quartimin(A)[[1]])
res$simplimax$f= f.simplimax(simplimax(A)[[1]])
res$geominQ$f= geominf(A_geo)
res$geomin0$f= geominf(Aprime,delta=0) 
res$promax$f=NA

mat=matrix(0,8,4)
colnames(mat)=c('mse','L0(c=0.01)','f(A)','f(rot(A))')
rownames(mat)=c('L1','L0.5','oblimin','quartimin','simplimax','geominQ','geomin(delat=0)','promax')
for (j in 1:8){
  for (m in 1:4){
    mat[j,m]=res[[j]][[m]]
  }
}

#setwd(path)
Alist = cbind(A,Aprime,A_geo,A_L1)
write.table(mat,file='mat1109_3.csv')
write.table(Alist,file='geomin_counterex1109_final.csv')





