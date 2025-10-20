# Simulation stuy for the paper "Parametric Modal Regression with Error Contaminated Covariates"
# Testing hypothesis: 
# Y~Gamma given link function and linear part are correctely specified.
#
# Test methods: MCCM1    

 rm(list=ls())
 packages = c("MASS", "smoothmest","matrixcalc","readxl","gofgamma")
 
 package.check = lapply( packages, FUN = function(x) 
    {
     if (!require(x, character.only = TRUE)) 
       {
         install.packages(x, dependencies = TRUE)
         library(x, character.only = TRUE)
       }
    } )

set.seed(6666)
# mode-link function
 
 g=function(x){exp(x)}

# Monte-Carlo Corrected log-likelihood function (mccl)
 
 mccl=function(par,mydata)
  {
   y=mydata$y
   data=mydata$data
   ZmB=cbind(data$ZmB.1,data$ZmB.2)
   impart=cbind(data$impart.Tp1,data$impart.Tp2)
   rp=ZmB%*%par[2:(d+1)]+par[1]
   ip=impart%*%par[2:(d+1)]
   Wb=complex(real=rp,imaginary=ip)
   logg=Re(apply(matrix(log(g(Wb)),nrow=B),2,mean))
   invg=Re(apply(matrix(1/g(Wb),nrow=B),2,mean))
   phi=par[d+2]
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
   rp=ZmB%*%par[2:(d+1)]+par[1]
   ip=impart%*%par[2:(d+1)]
   Wb=complex(real=rp,imaginary=ip)

   phi=par[d+2]

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
 Z1=data[,2:4]
 Z2=data[,5:7]
 y=scale(y,center=F,scale=T)
 Z1=scale(Z1,center=T,scale=T)
 Z2=scale(Z2,center=T,scale=T)

 Zrep=cbind(c(t(Z1)),c(t(Z2)))

 n=length(y)   # sample size
 B=50   # B-size
 nj=ncol(Z1)    # replication at each (X1,X2)
 d=2
   Zmean=matrix(rep(0,n*d),nrow=n)  # container of mean of Z's at each X1
   Zcov=matrix(rep(0,d*d*n),nrow=n,ncol=d*d) # container of covariance of Z's at each X1
   Zsd=matrix(rep(0,d*d*n),nrow=n,ncol=d*d) # container of covariance of Z's at each X1
 
   for(j in seq(n))
    {
     Zmean[j,]=apply(Zrep[(nj*(j-1)+1):(nj*j),],2,mean)
     Zcov[j,]=c(cov(Zrep[(nj*(j-1)+1):(nj*j),]))

     # Spectral decomposition of covariance
     Temp1=eigen(cov(Zrep[(nj*(j-1)+1):(nj*j),]))
     Temp2=(Temp1$vectors)%*%diag(sqrt(Temp1$values))%*%t((Temp1$vectors))
     Zsd[j,]=c(Temp2)
    }
 
   ZmB=Zmean[rep(seq_len(nrow(Zmean)),each=B),] # repeat each row of Zmean B times  
   ZvB=Zsd[rep(seq_len(nrow(Zsd)),each =B), ] # repeat each row of Zcov B times 
 
   Tpm=mvrnorm(n*B,rep(0,nj-1),diag(nj-1))      # Generating Tp-values 
   Tpf=function(vv){vv[1:d]/sqrt(sum(vv^2))}
   Tp=t(apply(Tpm,1,Tpf))        
   
   Tp1=apply(ZvB[,1:d]*Tp,1,sum)          # Compute S*Tp
   Tp2=apply(ZvB[,(d+1):(2*d)]*Tp,1,sum)
 
   impart=sqrt((nj-1)/nj)*cbind(Tp1,Tp2)  # compute the imaginary part 
   realpart=ZmB                           # extract the real part
 
 
   mydata=list(y=y,data=data.frame(ZmB=ZmB,impart=impart))
   bett=optim(par=c(-1.58,0.27,0.27,5),fn=mccl,mydata=mydata)$par
     
  Tv=mccl_T(par=bett,mydata=mydata)
  
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
  estt=matrix(0,nrow=M,ncol=4) 
  pvalue=rep(0,M)
  mm=0
  
  while(mm<M)
   {
    # Generating reproduce measurement error data tildZ
    
    tildZm=matrix(rep(0,n*d),n,d)     ### container of mean of tildZ's at each hatX1
    tildZsd=matrix(rep(0,n*d*d),n,d*d) ### container of covariance of tildZ's at each hatX1
    
    for(i in 1:n)
     {
      Si=matrix(Zcov[i,],nrow=2)
      tildu=mvrnorm(nj,rep(0,d),Si)         ##produce measurement error tildu
      reptildZ=t(matrix(rep(hatX1[i,],nj),2,nj))+tildu # repeat each row of hatX1 nj times
      tildZm[i,]=apply(reptildZ,2,mean)
      
      TildZ1=eigen(cov(reptildZ))
      TildZ2=(TildZ1$vectors)%*%diag(sqrt(TildZ1$values))%*%t((TildZ1$vectors))
      tildZsd[i,]=c(TildZ2)
    }
    
    # generating sample from the reproduce gamma mode regression 
    
    shapet=1+bett[d+2]
    ivegt=1/g(bett[1]+hatX1%*%bett[2:(d+1)])
    scalet=bett[d+2]*ivegt
    tildy=rgamma(n,shape=shapet,rate=scalet)
    
    ZmBt=tildZm[rep(seq_len(nrow(tildZm)),each=B),] # repeat each row of tildZm B times  
    ZvBt=tildZsd[rep(seq_len(nrow(tildZsd)),each =B), ] # repeat each row of tildZcov B times 
    
    Tpmt=mvrnorm(n*B,rep(0,nj-1),diag(nj-1))      # Generating Tp values 
    Tpft=function(vv){vv[1:d]/sqrt(sum(vv^2))}
    Tpt=t(apply(Tpmt,1,Tpft))        
    
    Tp1t=apply(ZvBt[,1:d]*Tpt,1,sum)          # Compute S*Tp
    Tp2t=apply(ZvBt[,(d+1):(2*d)]*Tpt,1,sum)
    
    impartt=sqrt((nj-1)/nj)*cbind(Tp1t,Tp2t)  # compute the imaginary part 
    realpartt=ZmBt                           # extract the real part
    
    mydata=list(y=tildy,data=data.frame(ZmB=ZmBt,impart=impartt))
    est=optim(par=c(-1.58,0.27,0.27,5),fn=mccl,mydata=mydata)
    
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




