# set path to be the current folder
setwd(paste0(path,'/Rdata'))
n_row=30
n_col=5
sim=500
x <-seq(0.01,0.5,by=0.01)
mse=matrix(0,length(x),4)
colnames(mse)<-paste0('N=',400*2^(0:3))
for ( j in 1:4){
  N<-400*2^(j-1)
  fn = paste0("ushape1012.unfixB.J",n_row,"K",n_col,"nsim", sim,"nsample", N, ".RData")
  load(fn)
  mse[,j] = rowMeans(sapply(1:sim,function(i,out){out[[i]]$mse},out=out),na.rm=T)
}
l1mse30=c(0.009235, 0.004462, 0.00225, 0.001137336,0.0005630748)
mse0=rbind(l1mse30[-5],mse)
pdf(file='unfixB1115_30S.pdf',width=8, height = 6)
matplot(c(0,x)[1:21],mse0[1:21,], type="l", lwd=2, xlab="Penalty Parameter", ylab= "MSE(Loadings)", axes=F)
axis(2)
axis(1)
points(rep(0,4),mse0[1,],type='p',col=1:4,pch=16)
legend(legend = colnames(mse),"topleft",col=1:4, lty=1:4, lwd=2)
dev.off()


n_row=15
n_col=3
sim=500
x <-seq(0.01,0.5,by=0.01)
mse=matrix(0,length(x),4)
colnames(mse)<-paste0('N=',400*2^(0:3))
for ( j in 1:4){
  N<-400*2^(j-1)
  fn = paste0("ushape1012.unfixB.J",n_row,"K",n_col,"nsim", sim,"nsample", N, ".RData")
  load(fn)
  mse[,j] = rowMeans(sapply(1:sim,function(i,out){out[[i]]$mse},out=out),na.rm=T)
}
l1mse15=c(0.010286, 0.005094, 0.002589,0.001300144,0.0006549897)
mse0=rbind(l1mse15[-5],mse)
pdf(file='unfixB1115_15S.pdf',width=8, height = 6)
matplot(c(0,x)[1:21],mse0[1:21,], type="l", lwd=2, xlab="Penalty Parameter", ylab= "MSE(Loadings)", axes=F)
axis(2)
axis(1)
points(rep(0,4),mse0[1,],type='p',col=1:4,pch=16)
legend(legend = colnames(mse),"topleft",col=1:4, lty=1:4, lwd=2)
dev.off()