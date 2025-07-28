# Simulation stuy for the paper "Parametric Modal Regression with Error Contaminated Covariates"
#install.packages("matrixcalc")
rm(list=ls())
ts=Sys.time()

# Load packages
library(MASS)
library(smoothmest)
library(pracma) #complex number
library(matrixcalc)

# mode-link function
g=function(x){exp(x)}

# Monte-Carlo Corrected log-likelihood function (mccl)

mccl=function(par,mydata)
 {
  y=mydata$y
  data=mydata$data
  ZmB=cbind(data$ZmB.1,data$ZmB.2)
  impart=cbind(data$impart.Tp1,data$impart.Tp2)
  X2rep=data$X2rep
  rp=ZmB%*%par[2:(d+1)]+as.vector(par[d+2]*X2rep)+par[1]
  ip=impart%*%par[2:(d+1)]
  Wb=complex(real=rp,imaginary=ip)
  logg=Re(apply(matrix(log(g(Wb)),nrow=B),2,mean))
  invg=Re(apply(matrix(1/g(Wb),nrow=B),2,mean))
  phi=par[d+3]
  out1=n*(1+phi)*log(phi)-n*log(gamma(phi+1))+phi*sum(log(y))
  out2=-(1+phi)*sum(logg)-phi*sum(y*invg)
  MCCL=out1+out2
  return(-MCCL/n)
}


# Function Test statistic T

mccl_T=function(par,mydata)
 {
  y=mydata$y
  data=mydata$data
  ZmB=cbind(data$ZmB.1,data$ZmB.2)
  impart=cbind(data$impart.Tp1,data$impart.Tp2)
  X2rep=data$X2rep
  phi=par[d+3]
  rp=ZmB%*%par[2:(d+1)]+as.vector(par[d+2]*X2rep)+par[1]
  ip=impart%*%par[2:(d+1)]
  Wb=complex(real=rp,imaginary=ip)
  
   logg=Re(apply(matrix(log(g(Wb)),nrow=B),2,mean))
   loginvg=Re(apply(matrix((g(Wb)/phi)*(log(phi/g(Wb))),nrow=B),2,mean))
   gphi=Re(apply(matrix(g(Wb)/phi,nrow=B),2,mean))

   Elogy=log(y)-digamma(1+phi)+log(phi)-logg
   Eylogy=y*log(y)-(1+phi)*digamma(2+phi)*gphi+(1+phi)*loginvg
  
  hatQ=cbind(Elogy,Eylogy)
  barQ=apply(hatQ,2,mean)
  hatsigma=solve(cov(hatQ)/n)
  return(t(barQ)%*%hatsigma%*%barQ)
}



set.seed(2222)
total=100
result=NULL 
Pvalue=NULL
jj=0
bett=matrix(rep(0,5*total),nrow=total)
while(jj<total)
{
  # Preparing Data
  n=100   # sample size
  B=100   # B-size
  d=2     # X1,X2 are contaminated with normal error
  nj=4    # replication at each (X1,X2)
  mu=rep(0,d)   # mean of measurement error U
  su2=1      # variance of measurement error U
  covx1=diag(d)  # covariance of (X1,X2)
  covu=diag(rep(su2,d))  # covariance of u
  X2<-rbinom(n,1,0.5)    # generate  n sample from X3
  X1=matrix(0,nrow=n,ncol=2)# generate n sample from (X1,X2)
  for(i in seq(n))
   {
    X1[i,]<-mvrnorm(1,rep((X2[i]==1)-(X2[i]==0),d),covx1)
   }

 bt0=-0.25                # intercept
 bt1=c(-0.25,-0.25)       # true beta-values for X1, X2
 bt2=-0.25                # true beta-value for X3
 phi=5                    # true phi-value
 bt3=-0.08

 # generating sample from the gamma mode regression 
 
 alph=1+phi               
 bt=phi/g(bt0+X1%*%bt1+X2*bt2+X1[,1]*X1[,1]*bt3)
 y=rgamma(n,shape=alph,rate=bt)

  # Generating normal measurement error
  
  U=mvrnorm(nj*n,mu,covu)  
  X1rep=X1[rep(seq(nrow(X1)),each =nj), ]  # repeat each row of X1 nj times
  Zrep=X1rep+U   # nj replicated observations at each row of X1
  
  Zmean=matrix(rep(0,n*d),nrow=n)  # container of mean of Z's at each X1
  Zcov=matrix(rep(0,d*d*n),nrow=n,ncol=d*d) # container of covariance of Z's at each X1
  Zsd=matrix(rep(0,d*d*n),nrow=n,ncol=d*d) # container of covariance of Z's at each X1
   
  for(j in seq(n))
   {
    Zmean[j,]=apply(Zrep[(nj*(j-1)+1):(nj*j),],2,mean)  #barZj
    Zcov[j,]=c(cov(Zrep[(nj*(j-1)+1):(nj*j),]))
    # Spectral decomposition of covariance
    Temp1=eigen(cov(Zrep[(nj*(j-1)+1):(nj*j),]))
    Temp2=(Temp1$vectors)%*%diag(sqrt(Temp1$values))%*%t((Temp1$vectors))  
    Zsd[j,]=c(Temp2)                                          ##Sj
   }
  
  ZmB=Zmean[rep(seq(n),each=B),] # repeat each row of Zmean B times  
  ZvB=Zsd[rep(seq(n),each =B),] # repeat each row of Zcov B times 
  
  Tpm=mvrnorm(n*B,rep(0,nj-1),diag(nj-1))      # Generating Td-values 
  Tpf=function(vv){vv[1:d]/sqrt(sum(vv^2))}
  Tp=t(apply(Tpm,1,Tpf))        
  
  Tp1=apply(ZvB[,1:d]*Tp,1,sum)          # Compute S*Td
  Tp2=apply(ZvB[,(d+1):(2*d)]*Tp,1,sum)
  
  impart=sqrt((nj-1)/nj)*cbind(Tp1,Tp2)  # compute the imaginary part 

  X2rep=rep(X2,each=B)
  

  mydata=list(y=y,data=data.frame(ZmB=ZmB,impart=impart,X2rep=X2rep))
  
  betts=optim(par=c(-0.3,-0.2,-0.2,-0.3,4),fn=mccl,mydata=mydata)
  

  
  if(betts$convergence==0)
   {
    point_estts = betts$par
    jj=jj+1
   } else 
   {
    next
   }
  
  bett=point_estts
  
  Tv=mccl_T(par=point_estts,mydata=mydata)
  
  # Bootstrapping Test Procedure
  
  # Generating hatX1
  
  barZm=apply(Zmean,2,mean)      ##sample mean of barZ
  barZmrep=matrix(rep(barZm,each=n),n,d) ##each row of barZ n times   
  biasZmbarZ=Zmean-barZmrep     ## sample bias  barZi-barZ
  hatsigmaz=cov(Zmean)          ###hatsigmaz
  
  # Spectral decomposition of covariance hatsigmaz
  Tempz1=eigen(solve(hatsigmaz))
  Tempz2=(Tempz1$vectors)%*%diag(sqrt(Tempz1$values))%*%t((Tempz1$vectors))
  hatsigmazz=Tempz2
  
  hatsigmau1=apply(Zcov,2,mean)
  hatsigmau=matrix(hatsigmau1,nrow=2)/nj   ##hatsigmau:barSi^2
  
  # choice sigmahat:(hatsigmaz-hatsigmau)+#
  hatsigmazu=round(hatsigmaz-hatsigmau,8)
  choicez=is.positive.definite(hatsigmazu)                 
  choiceu=is.positive.definite(-hatsigmazu)
  chatsigma=choicez*hatsigmazu+choiceu*0
  # Spectral decomposition of covariance chatsigma
  Tempc1=eigen(chatsigma)
  Tempc2=(Tempc1$vectors)%*%diag(sqrt(Tempc1$values))%*%t((Tempc1$vectors))
  choicsigmax=Tempc2
  
  ##sigmax/sigmaz
  choicsigma=choicsigmax%*%hatsigmazz
  
  # hatX1
  cz1=choicsigma[1,1]*biasZmbarZ[,1]+choicsigma[1,2]*biasZmbarZ[,2]
  cz2=choicsigma[2,1]*biasZmbarZ[,1]+choicsigma[2,2]*biasZmbarZ[,2]  
  hatX1=barZmrep+cbind(cz1,cz2)
  
  M=100
  estt=matrix(0,nrow=M,ncol=5) 
  pvalue=rep(0,M)
  mm=0
  
  while(mm<M)
   {
    # Generating reproduce measurement error data tiaotaZ
    
    tildZm=matrix(rep(0,n*d),n,d)     ### container of mean of tiaotaZ's at each hatX1
    tildZsd=matrix(rep(0,n*d*d),n,d*d) ### container of covariance of tiaotaZ's at each hatX1
    
    for(i in 1:n)
     {
      Si=matrix(Zcov[i,],nrow=2)
      tildu=mvrnorm(nj,rep(0,d),Si)         ##produce measurement error tiaotau
      reptildZ=t(matrix(rep(hatX1[i,],nj),2,nj))+tildu # repeat each row of hatX1 nj times
      tildZm[i,]=apply(reptildZ,2,mean)
      
      TildZ1=eigen(cov(reptildZ))
      TildZ2=(TildZ1$vectors)%*%diag(sqrt(TildZ1$values))%*%t((TildZ1$vectors))
      tildZsd[i,]=c(TildZ2)
    }
    
    # generating sample from the reproduce gamma mode regression 
    
    shapet=1+bett[d+3]
    ivegt=1/g(bett[1]+hatX1%*%bett[2:(d+1)]+X2*bett[d+2])
    scalet=bett[d+3]*ivegt
    tildy=rgamma(n,shape=shapet,rate=scalet)
    
    ZmBt=tildZm[rep(seq_len(nrow(tildZm)),each=B),] # repeat each row of tiaotaZm B times  
    ZvBt=tildZsd[rep(seq_len(nrow(tildZsd)),each =B), ] # repeat each row of tiaotaZcov B times 
    
    Tpmt=mvrnorm(n*B,rep(0,nj-1),diag(nj-1))      # Generating Tp values 
    Tpft=function(vv){vv[1:d]/sqrt(sum(vv^2))}
    Tpt=t(apply(Tpmt,1,Tpft))        
    
    Tp1t=apply(ZvBt[,1:d]*Tpt,1,sum)          # Compute S*Tp
    Tp2t=apply(ZvBt[,(d+1):(2*d)]*Tpt,1,sum)
    
    impartt=sqrt((nj-1)/nj)*cbind(Tp1t,Tp2t)  # compute the imaginary part 
    realpartt=ZmBt                           # extract the real part
    
    mydata=list(y=tildy,data=data.frame(ZmB=ZmBt,impart=impartt,X2rep=X2rep))
    est=optim(par=c(-0.3,-0.2,-0.2,-0.3,4),fn=mccl,mydata=mydata)
    
    if(est$convergence==0){
      point_estt = est$par
      mm=mm+1
    } else {
      next
    }
    
    estt[mm,]=point_estt
   
    Tvm=mccl_T(par=estt[mm,],mydata=mydata)
    
    
    pvalue[mm]=1*(Tvm>Tv)
    
    cat(mm,"\n")
  }
  estt
  pvalue
  
  mean(pvalue)
  result=rbind(result,estt)
  Pvalue=rbind(Pvalue,mean(pvalue))
  print(Pvalue)
  
  cat(jj,"\n")
}


mean(Pvalue<0.05)
mean(Pvalue<=0.05)

 te=Sys.time()
time=te-ts
save(file="n=100.misM2.RData",list=ls())


