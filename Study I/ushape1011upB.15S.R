rm(list=ls());
library(doParallel)
library(lavaan)
library(gtools)
library(zoo)
library(GPArotation)

setwd(path)
source('cl1121_Lp_fun.R')
sim = 500

lambda_list<-seq(0.01,0.5,by=0.01)
t1=proc.time()

newlasso.path<-function(lambda_list,L_vari,B,Psi_cfa,S,N,L){
  nlambda<-length(lambda_list)
  TPR = TNR = mse = TR = t = it = rep(0,nlambda)
  L_old = L != 0
  L0_list = list()
  t = it = 0
  L0 = L_vari
  B0 = diag(rep(1,n_col))
  Psi0 = log(Psi_cfa)
  for (j in 1:nlambda){
    ans_est <- prox_grad(L0,B0,Psi0,S,lambda_list[j],0)
    L0_list[[j]] = L0=permfilp(ans_est$L,L)
    Psi0 = ans_est$Psi
    B0 = ans_est$B
    t[j] = ans_est$t
    it[j] = ans_est$it
    #selection accuracy
    ans_hard=hard.each(0,L0,L_old)
    TPR[j] = ans_hard$TPR
    TNR[j] = ans_hard$TNR
    TR[j] = ans_hard$TR
    mse[j] = mean((L0-L)^2)
  }
    id <- 1:(length(TNR)+1)
    AUC <- sum(diff((1-c(1,TPR))[id])*rollmean(c(0,TNR)[id],2))
  return(list(L=L0_list,t=t,it=it, 
              TPR=TPR,TNR=TNR,TR=TR,AUC=AUC,
              mse=mse))
}

#15*3small loading##########################################
L<-matrix( c(0.71,0.00,0.00,
             0.00,0.75,0.00,
             0.00,0.00,0.83,
             0.96,0.00,0.00,
             0.00,0.68,0.00,
             0.00,0.00,0.96,
             0.98,0.00,0.00,
             0.00,0.86,0.00,
             0.00,0.00,0.85,
             0.62,0.35,0.00,
             0.00,0.68,0.42,
             0.50,0.00,0.67,
             0.87,0.00,0.31,
             0.43,0.75,0.00,
             0.00,0.48,0.91),nrow=15,ncol=3,byrow<-TRUE);
B<-matrix(c(1,0.02079577,0.5024378,
            0,0.99978374,0.2635086,
            0,0,0.8234801),3,3,byrow=TRUE)
P<-c(0.24, 0.32, 0.45, 0.65, 0.18, 0.64, 0.67, 0.51, 0.49, 0.06, 0.19, 0.16, 0.52, 0.33, 0.57)
##############################################################################
n_row<- nrow(L)
n_col<- ncol(L)
Sig <- L%*%t(B)%*%B%*%t(L)+diag(exp(P))
(Sig)
(exp(P))
(t(B)%*%B)
L_old=L!=0
seed=0

############ Replicates
ncores=detectCores()
cl = makeCluster(ncores-1)
registerDoParallel(cl)


for ( j in 1:5){
  t=proc.time()
  N<-400*2^(j-1)
out = foreach(i = 1:sim, .packages = c("lavaan","gtools","zoo","GPArotation"), .errorhandling = "pass") %dopar%{
  set.seed(i)
  res=list()
  ## generate data
  SLP0=Gen.SLP(Sig,N,n_col)
  S=SLP0$S
  L_vari=SLP0$L_vari
  Psi_cfa=SLP0$Psi_cfa
  
  ans_lasso<-newlasso.path(lambda_list,L_vari,B,Psi_cfa,S,N,L)
  res$L = ans_lasso$L
  res$t = ans_lasso$t
  res$it = ans_lasso$it
  res$mse = ans_lasso$mse
  #hard-thresholding
  res$TPR = ans_lasso$TPR
  res$TNR = ans_lasso$TNR
  res$TR = ans_lasso$TR
  res$AUC = ans_lasso$AUC
  return(res)
}
save(out,L,B,P,sim,N,seed, t,file = paste0("ushape1012.unfixB.J",nrow(L),"K",n_col,"nsim", sim,"nsample", N, ".RData"))
print(round((proc.time()-t)[3]/60,2))
}
stopCluster(cl)
