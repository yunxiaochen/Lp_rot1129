rm(list=ls());
library(doParallel)
library(lavaan)
library(gtools)
library(zoo)
library(numDeriv)
# 1114 
# only change is to use cl1110_Lp_fun.R
source('cl1121_Lp_fun.R')
sim = 500
N = 3000
c_list<-c(seq(0,0.19,by=0.01),seq(0.2,2,by=0.1),3)
c_bic_list<-c(seq(0.02,0.4,by=0.02))
lambda_list<-c(seq(0.01,0.37,by=0.01),0.4,0.5,2)
t1=proc.time()



#15*3small loading##########################################
# L<-matrix( c(0.5310173,0.99539848,0.0000000,
#              0.7442478,0.00000000,1.1991317,
#              0.0000000,1.98381219,0.9870826,
#              1.8164156,0.76007036,0.0000000,
#              0.4033639,0.00000000,1.6547466,
#              0.0000000,1.86941046,1.3369335,
#              1.8893505,0.42428504,0.0000000,
#              1.3215956,0.00000000,0.2158873,
#              0.0000000,0.25111019,1.4474219,
#              0.1235725,0.53444134,0.0000000,
#              0.4119491,0.00000000,1.6418926,
#              0.0000000,0.02678067,1.2941204,
#              1.3740457,0.76477591,0.0000000,
#              0.7682074,0.00000000,1.1060726,
#              0.0000000,0.68069799,1.0594392),nrow=15,ncol=3,byrow<-TRUE);
#> set.seed(1)
#L=matrix(2*runif(54),18,3)
L<-matrix(c(0.5310173, 0.76007036, 0.0000000,
            0.7442478,0.00000000,0.2158873
            ,0.0000000,1.86941046,1.4474219
            ,1.8164156,0.42428504,0.0000000
            ,0.4033639,0.00000000,1.6418926
            ,0.0000000,0.25111019,1.2941204
            ,1.8893505,0.53444134,0.0000000
            ,1.3215956,0.00000000,1.1060726
            ,0.0000000,0.02678067,1.0594392
            ,0.1235725,0.76477591,0.0000000
            ,0.4119491,0.00000000,0.0466624
            ,0.0000000,0.68069799,0.9544601
            ,1.3740457,0.96416023,0.0000000
            ,0.7682074,0.00000000,1.3854631
            ,0.0000000,0.98708261,0.9552392
            ,0.9953985,0.37243520,0.0000000
            ,1.4352370,0.00000000,0.8761942
            ,0.0000000,1.33693348,0.4895946),nrow=18,ncol=3,byrow<-TRUE)

# L<-matrix(c(1.558740486,0.85441467,0.0000000
# ,1.657020284,0.00000000,0.9317462
# ,0.000000000,0.58549927,1.2020446
# ,0.003862898,1.31511640,0.0000000
# ,1.325770219,0.00000000,1.3067797
# ,0.000000000,1.97972289,0.3759713
# ,1.747234310,1.80369291,0.0000000
# ,0.419183761,0.00000000,0.8012482
# ,0.000000000,1.42552833,0.8864383
# ,1.111660344,0.12507188,0.0000000
# ,1.690253911,0.00000000,0.6219600
# ,0.000000000,0.05416751,0.5556852
# ,1.747065776,1.30126546,0.0000000
# ,1.869674556,0.00000000,1.6040512
# ,0.000000000,1.79164432,0.1484092
# ,0.043441413,1.33081235,0.0000000
# ,1.615220623,0.00000000,0.8959223
# ,0.000000000,1.71232147,0.7156304
# ,0.264148903,0.94818884,0.0000000
# ,0.718479480,0.00000000,0.8437018
# ,0.000000000,1.24025397,0.6761058
# ,1.025781205,1.21810072,0.0000000
# ,0.200434512,0.00000000,1.5926942
# ,0.000000000,0.01061715,1.6727535
# ,0.898425022,1.63901289,0.0000000
# ,1.414164154,0.00000000,1.7692894
# ,0.000000000,1.42603030,0.7146371
# ,0.984604014,1.34992676,0.0000000
# ,0.766891322,0.00000000,1.6979906
# ,0.000000000,0.11207141,1.3625525
# ),nrow=30,ncol=3,byrow<-TRUE)
# #B<-diag(3)
B<-matrix(c(1,0.02079577,0.5024378,
            0,0.99978374,0.2635086,
            0,0,0.8234801),3,3,byrow=TRUE) #


####2 irls_ep=0.01 0.1,0.001--- 3 (1,1.5)(0.5,1) 4 sample size

P<-c(0.24, 0.32, 0.45, 0.65, 0.18, 0.64, 0.67, 0.51, 0.49, 0.06, 0.19, 0.16, 0.52, 0.33, 0.57,0.40,0.54,0.69)
#P<-c(0.24,0.32,0.45,0.65,0.18,0.64,0.67,0.51,0.49,0.06,0.19,0.16,0.52,0.33,0.57,0.40,0.54,0.69,0.32,0.58,0.66,0.19,0.50,0.12,0.24,0.33,0.01,0.32,0.63,0.29)

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

out = foreach(i = 1:sim, .packages = c("lavaan","gtools","zoo"), .errorhandling = "pass") %dopar%{
  set.seed(i)
  
  ## generate data
  SLP0=Gen.SLP(Sig,N,n_col)
  S=SLP0$S
  L_vari=SLP0$L_vari
  Psi_cfa=SLP0$Psi_cfa
  Psi0=log(Psi_cfa)
  
  res=list()
 # L_vari=permfilp(L_vari,L)
  #res$vari$mse= mean((L_vari-L)^2)
  
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
  return(res)
}
t=proc.time()[3]-t1[3]
save(out,L,B,P,sim,N,seed, t,file = paste0("1121Counter.18S",".RData"))

stopCluster(cl)