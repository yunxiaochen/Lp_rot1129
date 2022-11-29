rm(list=ls());
# changes made
# oblique rotations
# N=400,800,1600
# file-name
# replace the sgn(refit.bic,refit.L,CI) and CI function
library(doParallel)
library(lavaan)
library(gtools)
library(zoo)
library(GPArotation)
library(numDeriv)


source('cl1121_Lp_fun.R')
sim = 500

c_list<-c(seq(0,0.19,by=0.01),seq(0.2,2,by=0.1),3)
c_bic_list<-c(seq(0.02,0.4,by=0.02))
lambda_list<-c(seq(0.01,0.37,by=0.01),0.4,0.5,2)
t1=proc.time()


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
#########################################################
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



for (N in c(400,800,1600)){
out = foreach(i = 1:sim, .packages = c("lavaan","gtools","zoo","GPArotation","numDeriv"), .errorhandling = "pass") %dopar%{
  set.seed(i)
  
  ## generate data
  SLP0=Gen.SLP(Sig,N,n_col)
  S=SLP0$S
  L_vari=SLP0$L_vari
  Psi_cfa=SLP0$Psi_cfa
  Psi0=log(Psi_cfa)
  
  res=list()
  res$vari$L=L_vari #permfilp(L_vari,L)
  res$vari$mse= mean((oblimin(L_vari)[[1]]-L)^2)
  res$oblimin$mse= mean((oblimin(L_vari)[[1]]-L)^2)
  res$quartimin$mse= mean((quartimin(L_vari)[[1]]-L)^2)
  res$simplimax$mse= mean((simplimax(L_vari)[[1]]-L)^2)
  res$geominQ$mse= mean((geominQ(L_vari)[[1]]-L)^2)
  res$promax$mse=mean((promax(L_vari)[[1]]-L)^2)
  
  ## irls p=1
  ans_est1<-irls(1,L_vari)
  respefl=pefl.LT(ans_est1$L,L,ans_est1$T)
  L_est=respefl$L
  T_est=respefl$B

  res$irls1$L=L_est
  res$irls1$t=ans_est1$t
  res$irls1$it=ans_est1$it
  res$irls1$mse= mean((L_est-L)^2)
  #hard-thresholding
  ans_hard=hard(c_list,L_est,L_old)
  res$irls1$TPR =ans_hard$TPR
  res$irls1$TNR =ans_hard$TNR
  #res$irls1$AUC =ans_hard$AUC
  #bic-accuracy
  bic=sapply(c_bic_list,sgnrefit_optim2B2,L_irls=L_est,S=S,N=N,mod=0,k=-1,Psi0=Psi0,B=T_est)
  res$irls1$c= c0 = c_bic_list[which.min(bic)[1]]
  res$irls1$L_bic=L_bic=sgnrefit_optim2B2(c0,L_est,S,N,mod=1,Psi0=Psi0,B=T_est)$L
  res$irls1$L_bic.res=hard.each(c0,L_est,L_old)
  #CI-coverage
  ans_ci=CI3(c0,L_est,S,N,Psi0, T_est)
  res$irls1$L_class=ans_ci$L_class
  res$irls1$ci.accuracy=ans_ci$ci.accuracy
  res$irls1$L_lower=ans_ci$L_lower
  res$irls1$L_upper=ans_ci$L_upper
  res$irls1$na.flag=ans_ci$na.flag
  
  #irls p=0.5, the reason is initial value
  ans_est<-irls(0.5,L_vari,T=ans_est1$T)
  respefl=pefl.LT(ans_est$L,L,ans_est$T)
  L_est=respefl$L
  T_est=respefl$B
  
  res$irls0.5$L=L_est
  res$irls0.5$t=ans_est$t
  res$irls0.5$it=ans_est$it
  res$irls0.5$mse= mean((L_est-L)^2)
  #hard-thresholding
  ans_hard=hard(c_list,L_est,L_old)
  res$irls0.5$TPR =ans_hard$TPR
  res$irls0.5$TNR =ans_hard$TNR
  #res$irls0.5$AUC =ans_hard$AUC
  #bic-accuracy
  bic=sapply(c_bic_list,sgnrefit_optim2B2,L_irls=L_est,S=S,N=N,mod=0,k=-1,Psi0=Psi0,B=T_est)
  res$irls0.5$c= c0 = c_bic_list[which.min(bic)[1]]
  res$irls0.5$L_bic=L_bic=sgnrefit_optim2B2(c0,L_est,S,N,mod=1,Psi0=Psi0,B=T_est)$L
  res$irls0.5$L_bic.res=hard.each(c0,L_est,L_old)
  #CI-coverage
  ans_ci=CI3(c0,L_est,S,N,Psi0,T_est)
  res$irls0.5$L_class=ans_ci$L_class
  res$irls0.5$ci.accuracy=ans_ci$ci.accuracy
  res$irls0.5$L_lower=ans_ci$L_lower
  res$irls0.5$L_upper=ans_ci$L_upper
  res$irls0.5$na.flag=ans_ci$na.flag
  
  
  if(T){
    ## lasso.ap
    ans_est2<-prox_grad(L_vari,diag(rep(1,n_col)),log(Psi_cfa),S,0.01)
    respefl=pefl.LT(ans_est2$L,L,ans_est2$B)
    L_est=respefl$L
    T_est=respefl$B
    res$lasso.ap0.01$L=L_est
    res$lasso.ap0.01$t=ans_est2$t
    res$lasso.ap0.01$it=ans_est2$it
    res$lasso.ap0.01$mse= mean((L_est-L)^2)
    #hard-thresholding
    ans_hard=hard(c_list,L_est,L_old)
    res$lasso.ap0.01$TPR =ans_hard$TPR
    res$lasso.ap0.01$TNR =ans_hard$TNR
    #res$lasso.ap$AUC =ans_hard$AUC
    #bic-accuracy
    bic=sapply(c_bic_list,sgnrefit_optim2B2,L_irls=L_est,S=S,N=N,mod=0,k=-1,Psi0=ans_est2$Psi,B=rho(T_est+0.01*diag(n_col))) #,B=T_est
    res$lasso.ap0.01$c= c0 = c_bic_list[which.min(bic)[1]]
    res$lasso.ap0.01$L_bic=L_bic=sgnrefit_optim2B2(c0,L_est,S,N,mod=1,Psi0=ans_est2$Psi,B=rho(T_est+0.01*diag(n_col)))$L
    res$lasso.ap0.01$L_bic.res=hard.each(c0,L_est,L_old)
    
  
    ## lasso.ap
    ans_est3<-prox_grad(ans_est2$L,ans_est2$B,ans_est2$Psi,S,0.05)
    respefl=pefl.LT(ans_est3$L,L,ans_est3$B)
    L_est=respefl$L
    T_est=respefl$B
    res$lasso.ap0.05$L=L_est
    res$lasso.ap0.05$t=ans_est3$t
    res$lasso.ap0.05$it=ans_est3$it
    res$lasso.ap0.05$mse= mean((L_est-L)^2)
    #hard-thresholding
    ans_hard=hard(c_list,L_est,L_old)
    res$lasso.ap0.05$TPR =ans_hard$TPR
    res$lasso.ap0.05$TNR =ans_hard$TNR
    #res$lasso.ap0.05$AUC =ans_hard$AUC
    #bic-accuracy
    bic=sapply(c_bic_list,sgnrefit_optim2B2,L_irls=L_est,S=S,N=N,mod=0,k=-1,Psi0=ans_est3$Psi,B=rho(T_est+0.01*diag(n_col))) 
    res$lasso.ap0.05$c= c0 = c_bic_list[which.min(bic)[1]]
    res$lasso.ap0.05$L_bic=L_bic=sgnrefit_optim2B2(c0,L_est,S,N,mod=1,Psi0=ans_est3$Psi,B=rho(T_est+0.01*diag(n_col))) $L
    res$lasso.ap0.05$L_bic.res=hard.each(c0,L_est,L_old)
    
    
    ## lasso.ap
    ans_est4<-prox_grad(ans_est3$L,ans_est3$B,ans_est3$Psi,S,0.1)
    respefl=pefl.LT(ans_est4$L,L,ans_est4$B)
    L_est=respefl$L
    T_est=respefl$B
    res$lasso.ap0.1$L=L_est
    res$lasso.ap0.1$t=ans_est4$t
    res$lasso.ap0.1$it=ans_est4$it
    res$lasso.ap0.1$mse= mean((L_est-L)^2)
    #hard-thresholding
    ans_hard=hard(c_list,L_est,L_old)
    res$lasso.ap0.1$TPR =ans_hard$TPR
    res$lasso.ap0.1$TNR =ans_hard$TNR
    #res$lasso.ap0.1$AUC =ans_hard$AUC
    #bic-accuracy
    bic=sapply(c_bic_list,sgnrefit_optim2B2,L_irls=L_est,S=S,N=N,mod=0,k=-1,Psi0=ans_est4$Psi,B=rho(T_est+0.01*diag(n_col))) 
    res$lasso.ap0.1$c= c0 = c_bic_list[which.min(bic)[1]]
    res$lasso.ap0.1$L_bic=L_bic=sgnrefit_optim2B2(c0,L_est,S,N,mod=1,Psi0=ans_est4$Psi,B=rho(T_est+0.01*diag(n_col)))$L
    res$lasso.ap0.1$L_bic.res=hard.each(c0,L_est,L_old)
    
    ## lasso.ap
    ans_est5<-prox_grad(ans_est4$L,ans_est4$B,ans_est4$Psi,S,0.2)
    respefl=pefl.LT(ans_est5$L,L,ans_est5$B)
    L_est=respefl$L
    T_est=respefl$B
    res$lasso.ap0.2$L=L_est
    res$lasso.ap0.2$t=ans_est5$t
    res$lasso.ap0.2$it=ans_est5$it
    res$lasso.ap0.2$mse= mean((L_est-L)^2)
    #hard-thresholding
    ans_hard=hard(c_list,L_est,L_old)
    res$lasso.ap0.2$TPR =ans_hard$TPR
    res$lasso.ap0.2$TNR =ans_hard$TNR
    #res$lasso.ap0.2$AUC =ans_hard$AUC
    #bic-accuracy
    bic=sapply(c_bic_list,sgnrefit_optim2B2,L_irls=L_est,S=S,N=N,mod=0,k=-1,Psi0=ans_est5$Psi,B=rho(T_est+0.01*diag(n_col))) 
    res$lasso.ap0.2$c= c0 = c_bic_list[which.min(bic)[1]]
    res$lasso.ap0.2$L_bic=L_bic=sgnrefit_optim2B2(c0,L_est,S,N,mod=1,Psi0=ans_est5$Psi,B=rho(T_est+0.01*diag(n_col)))$L
    res$lasso.ap0.2$L_bic.res=hard.each(c0,L_est,L_old)
  
    ans_est6<-prox_grad(ans_est5$L,ans_est5$B,ans_est5$Psi,S,0.5)
    respefl=pefl.LT(ans_est6$L,L,ans_est6$B)
    L_est=respefl$L
    T_est=respefl$B
    res$lasso.ap0.5$L=L_est
    res$lasso.ap0.5$t=ans_est6$t
    res$lasso.ap0.5$it=ans_est6$it
    res$lasso.ap0.5$mse= mean((L_est-L)^2)
    #hard-thresholding
    ans_hard=hard(c_list,L_est,L_old)
    res$lasso.ap0.5$TPR =ans_hard$TPR
    res$lasso.ap0.5$TNR =ans_hard$TNR
    #res$lasso.ap0.5$AUC =ans_hard$AUC
    #bic-accuracy
    bic=sapply(c_bic_list,sgnrefit_optim2B2,L_irls=L_est,S=S,N=N,mod=0,k=-1,Psi0=ans_est6$Psi,B=rho(T_est+0.01*diag(n_col))) 
    res$lasso.ap0.5$c= c0 = c_bic_list[which.min(bic)[1]]
    res$lasso.ap0.5$L_bic=L_bic=sgnrefit_optim2B2(c0,L_est,S,N,mod=1,Psi0=ans_est6$Psi,B=rho(T_est+0.01*diag(n_col)))$L
    res$lasso.ap0.5$L_bic.res=hard.each(c0,L_est,L_old)
  }
  return(res)
}
t=proc.time()[3]-t1[3]
save(out,L,B,P,sim,N,seed, t, file = paste0("Lp.1121.",nrow(L),"N", N, ".RData"))
}
stopCluster(cl)