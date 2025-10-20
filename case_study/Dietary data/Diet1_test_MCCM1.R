# Simulation stuy for the paper "Parametric Modal Regression with Error Contaminated Covariates"
# Testing hypothesis: 
# Y~Gamma given link function and linear part are correctely specified.
#
# Test methods: MCCM1   
# Full data(D1), corresponding  Z=data[,2:7]  y=y
   #D2,  Z=data[y<5000,2:7]  y=y[y<5000]
   #D3,  Z=data[y<4000,2:7] y=y[y<4000]

rm(list=ls())
ts=Sys.time()

# Load packages
library(MASS)
library(smoothmest)
library(pracma) #complex number
library(matrixcalc)


set.seed(6666)

# mode-link function
g=function(x){exp(x)}

# Monte-Carlo Corrected log-likelihood function (mccl)

mccl=function(par,mydata)
 {
  y=mydata$y
  data=mydata$data
  ZmB=data$ZmB
  impart=data$impart.Tp
  rp=ZmB*par[2]+par[1]
  ip=impart*par[2]
  Wb=complex(real=rp,imaginary=ip)
  logg=Re(apply(matrix(log(g(Wb)),nrow=B),2,mean))
  invg=Re(apply(matrix(1/g(Wb),nrow=B),2,mean))
  phi=par[3]
  out1=n*(1+phi)*log(phi)-n*log(gamma(phi+1))+phi*sum(log(y))
  out2=-(1+phi)*sum(logg)-phi*sum(y*invg)
  MCCL=out1+out2
  return(-MCCL/n)
}


# Function Test statistic T

mccl_T=function(par,mydata)
 {
   phi=par[3]
   y=mydata$y
   data=mydata$data
   ZmB=data$ZmB
   impart=data$impart.Tp
   rp=ZmB*par[2]+par[1]
   ip=impart*par[2]
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


 data<- read.csv("wishreg.csv", header=T)

 y=data$ffq
 Z=data[y<5000,2:7]
 y=y[y<5000]
 y=scale(y,center=F,scale=T)
 Z=scale(Z,center=T,scale=T)
 Zrep=c(t(Z))

 n=length(y)   # sample size
 B=50   # B-size
 nj=ncol(Z)    # replication at each (X1,X2)



 Zmean=rep(0,n)  # container of mean of Z's at each X1
 Zsd=rep(0,n) # container of sd of Z's at each X1
 Zcov=rep(0,n) # container of cov of Z's at each X1

 for(j in seq(n))
   {
    Zmean[j]=mean(Zrep[(nj*(j-1)+1):(nj*j)])
    Zcov[j]=var(Zrep[(nj*(j-1)+1):(nj*j)])
    # Spectral decomposition of covariance
    Zsd[j]=sqrt(var(Zrep[(nj*(j-1)+1):(nj*j)]))
   }

 ZmB=rep(Zmean,each=B) # repeat each row of Zmean B times  
 ZvB=rep(Zsd,each =B)  # repeat each row of Zcov B times 
 
 Tpm=mvrnorm(n*B,rep(0,nj-1),diag(nj-1))      # Generating Tp-values 
 Tpf=function(vv){vv[1:1]/sqrt(sum(vv^2))}
 Tp=t(apply(Tpm,1,Tpf))        
   
 Tp=ZvB*Tp          # Compute S*Tp
 
 impart=sqrt((nj-1)/nj)*Tp   # compute the imaginary part 
 realpart=ZmB                           # extract the real part
 
 mydata=list(y=y,data=data.frame(ZmB=ZmB,impart=impart))
  
 betts=optim(par=c(-1.58,0.27,5),fn=mccl,mydata=mydata)

 Tv=mccl_T(par=betts$par,mydata=mydata)


  # Bootstrapping Test Procedure
  
  # Generating hatX1
  
  barZm=mean(Zmean)      ##sample mean of barZ
  barZmrep=rep(barZm,each=n) ##each row of barZ n times   
  biasZmbarZ=Zmean-barZmrep     ## sample bias  barZi-barZ
  hatsigmaz=var(Zmean)          ###hatsigmaz
  
  # Spectral decomposition of covariance hatsigmaz
   hatsigmazz=sqrt(solve(hatsigmaz))
  
  hatsigmau=mean(Zcov)/nj
  
  # choice sigmahat:(hatsigmaz-hatsigmau)+#
  hatsigmazu=round(hatsigmaz-hatsigmau,8)
  chatsigma=hatsigmazu*(hatsigmazu>0)+(hatsigmazu<0)*0
  # Spectral decomposition of covariance chatsigma
  choicsigmax=sqrt(chatsigma)
  
  ##sigmax/sigmaz
  choicsigma=choicsigmax*hatsigmazz
  
  # hatX1
  hatX1=barZmrep+c(choicsigma)*biasZmbarZ

  M=100
  estt=matrix(0,nrow=M,ncol=3) 
  pvalue=rep(0,M)
  mm=0
  
  while(mm<M)
   {
    # Generating reproduce measurement error data tildZ
    
    tildZm=rep(0,n)     ### container of mean of tildZ's at each hatX1
    tildZsd=rep(0,n) ### container of covariance of tildZ's at each hatX1
    
    for(i in 1:n)
     {
      Si=Zcov[i]
      tildu=rnorm(nj,0,Si)         ##produce measurement error tildu
      reptildZ=rep(hatX1[i],each=nj)+tildu # repeat each row of hatX1 nj times
      tildZm[i]=mean(reptildZ)
      tildZsd[]=sqrt(var(reptildZ))
    }
    
    # generating sample from the reproduce gamma mode regression 
    bett=betts$par
    shapet=1+bett[3]
    ivegt=1/g(bett[1]+hatX1*bett[2])
    scalet=bett[3]*ivegt
    tildy=rgamma(n,shape=shapet,rate=scalet)
    
    ZmBt=rep(tildZm,each=B) # repeat each row of tildZm B times  
    ZvBt=rep(tildZsd,each =B) # repeat each row of tildZcov B times 
    
    Tpmt=mvrnorm(n*B,rep(0,nj-1),diag(nj-1))      # Generating Tp values 
    Tpft=function(vv){vv[1:1]/sqrt(sum(vv^2))}
    Tpt=t(apply(Tpmt,1,Tpft))        
    
    Tpt=ZvBt*Tpt          # Compute S*Tp
    
    impartt=sqrt((nj-1)/nj)*Tpt  # compute the imaginary part 
    realpartt=ZmBt                           # extract the real part
    
    mydata=list(y=tildy,data=data.frame(ZmB=ZmBt,impart=impartt))
    est=optim(par=c(-1.58,0.27,5),fn=mccl,mydata=mydata)
    
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
  


