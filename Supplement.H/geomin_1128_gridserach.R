# this part contains the grid search for L1 solution, the calculation requires a long time.
gridlen=0.05

# Generate A
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

# grid serach of the maximum value of L1,L0.5 over the whole oblique rotation matrix
n <- 9
l <- rep(list(0:1), n)
sign=as.matrix(2*expand.grid(l)-1)


i=0
for( t11 in seq(0,1,by = gridlen)){
  for (t21 in seq(0,sqrt(1-t11^2),by=gridlen)){
    i=i+1
  }
}
length=i

# To use the parallel computing, we devide the task of 
# going through the whole oblique rot matrix space by different first column entries.
tlist=matrix(0,length,2)
i=0
for( t11 in seq(0,1,by = gridlen)){
  for (t21 in seq(0,sqrt(1-t11^2),by=gridlen)){
    tlist[i,]=c(t11,t21)
    i=i+1
  }
}

ncores=detectCores()
cl = makeCluster(ncores-1)
registerDoParallel(cl)
for (p in c(1,0.5)){
  t1=proc.time()[3]
  out = foreach(i = 1:nrow(tlist),.packages='matrixcalc', .errorhandling = "pass") %dopar%{
    minQ=Inf
    res=list()
    t11=tlist[i,1]
    t21=tlist[i,2]
    t31=sqrt(1-t11^2-t21^2)
    for( t12 in seq(0,1,by = gridlen)){
      for (t22 in seq(0,sqrt(1-t12^2),by=gridlen)){
        t32=sqrt(1-t12^2-t22^2)
        for( t13 in seq(0,1,by = gridlen)){
          for (t23 in seq(0,sqrt(1-t13^2),by=gridlen)){
            t33=sqrt(1-t13^2-t23^2)
            for(j in 1:nrow(sign)){
              Tmat=matrix(sign[j,]*c(t11,t21,t31,t12,t22,t32,t13,t23,t33),3,3)
              if (!is.singular.matrix(Tmat)){
                valQ=sum(abs(A%*%solve(t(Tmat)))^p)
                if(valQ<minQ){
                  minQ=valQ
                  minT=Tmat
                }}
            }
          }}
      }}
    res$minQ=minQ
    res$minT=minT
    return(res)
  }
  print(proc.time()[3]-t1)
  save(out,file = paste0("1128A_L",p,".RData"))
}
stopCluster(cl)
