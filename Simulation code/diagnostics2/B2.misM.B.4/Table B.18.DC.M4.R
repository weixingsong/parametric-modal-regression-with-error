# Simulation stuy for the paper "Parametric Modal Regression with Error Contaminated Covariates"
#install.packages("matrixcalc")
rm(list=ls())
ts=Sys.time()

# Load packages
 library(MASS)
 library(smoothmest)
 library(pracma) #complex number
 library(dcov)   #correlation of distances
 library(matrixcalc)

# mode-link function
 g=function(x){1/(1+exp(x))}

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

 Impart=function(n,ZvB)
  {
   Tpm=mvrnorm(n*B,rep(0,nj-1),diag(nj-1))      # Generating Td-values 
   Tpf=function(vv){vv[1:d]/sqrt(sum(vv^2))}
   Tp=t(apply(Tpm,1,Tpf))        
  
   Tp1=apply(ZvB[,1:d]*Tp,1,sum)          # Compute S*Td
   Tp2=apply(ZvB[,(d+1):(2*d)]*Tp,1,sum)
  
   impart=sqrt((nj-1)/nj)*cbind(Tp1,Tp2)  # compute the imaginary part 
   return(impart)
  }


 RTn=function(par,mydata1,mydata2)
  {
   y1=mydata1$y
   y2=mydata2$y

   data1=mydata1$data
   data2=mydata2$data

   ZmB1=cbind(data1$ZmB.1,data1$ZmB.2)
   ZmB2=cbind(data2$ZmB.1,data2$ZmB.2)

   impart1=cbind(data1$impart.Tp1,data1$impart.Tp2)
   impart2=cbind(data2$impart.Tp1,data2$impart.Tp2)

   X2rep1=data1$X2rep
   X2rep2=data2$X2rep

   hatphi=par[d+3]

   rp1=ZmB1%*%par[2:(d+1)]+as.vector(par[d+2]*X2rep1)+par[1]
   rp2=ZmB2%*%par[2:(d+1)]+as.vector(par[d+2]*X2rep2)+par[1]

   ip1=impart1%*%par[2:(d+1)]
   ip2=impart2%*%par[2:(d+1)]

   Wb1=complex(real=rp1,imaginary=ip1)
   Wb2=complex(real=rp2,imaginary=ip2)

   invg1=Re(apply(matrix(1/g(Wb1),nrow=B),2,mean))
   invg2=Re(apply(matrix(1/g(Wb2),nrow=B),2,mean))
   invgivg=Re(apply(matrix(g(Wb2),nrow=B),2,mean))

   R1=hatphi*y1*invg1+hatphi*y2*invg2
   R2=(y1/y2)*invg1*invgivg

 ##check the independence between R1 and R2
  S21<-(4/n^2)*sum(abs(R1%x%rep(1,n/2)-rep(R1,n/2))) 
  S22<-(4/n^2)*sum(abs(R2%x%rep(1,n/2)-rep(R2,n/2)))    #compute S2
  S2=S21*S22
  Tn=dcov(R1,R2)
  return(sqrt((n*Tn)/(2*S2)))
 }

set.seed(555)
total=100
Pvalue=NULL
jj=0

while(jj<total)
{
  # Preparing Data
  n=200   # sample size
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
  m=5                    # true phi-value

  # generating sample from the gamma mode regression 
 
 thetax=bt0+X1%*%bt1+X2*bt2
 alph=1+m*g(thetax)
 bt=1+m*(1-g(thetax))
 y=rbeta(n,shape1=alph,shape2=bt)


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
  
  ZmB1=ZmB[1:(n*B/2),]
  ZmB2=ZmB[(n*B/2+1):(n*B),]

  ZvB1=ZvB[1:(n*B/2),]
  ZvB2=ZvB[(n*B/2+1):(n*B),]

  impart=Impart(n,ZvB)
  impart1=Impart(n/2,ZvB1)
  impart2=Impart(n/2,ZvB2)

  X2rep=rep(X2,each=B)       # repeat each row of X2 B times
  X2rep1=rep(X2[1:(n/2)],each=B)
  X2rep2=rep(X2[(n/2+1):n],each=B)

  y1=y[1:(n/2)]
  y2=y[(n/2+1):n]
  

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

  data1=list(y=y1,data=data.frame(ZmB=ZmB1,impart=impart1,X2rep=X2rep1))
  data2=list(y=y2,data=data.frame(ZmB=ZmB2,impart=impart2,X2rep=X2rep2))

  
  Tv=RTn(par=point_estts,mydata1=data1,mydata2=data2)
  
 # hatX1
 barZm=apply(Zmean,2,mean)      ##sample mean of barZ
 barZmrep=t(matrix(rep(barZm,n),d,n)) ##each row of barZ n times   
 biasZmbarZ=Zmean-barZmrep     ## sample bias  barZi-barZ
 hatsigmaz=cov(biasZmbarZ)          ###hatsigmaz

 # Spectral decomposition of covariance hatsigmaz
 Tempz1=eigen(solve(hatsigmaz))
 Tempz2=(Tempz1$vectors)%*%diag(sqrt(Tempz1$values))%*%t((Tempz1$vectors))
 hatsigmazz=Tempz2

 hatsigmau1=apply(Zcov,2,mean)
 hatsigmau=matrix(hatsigmau1,nrow=2)/nj   ##hatsigmau:barSi^2

 # choice sigmahat:(hatsigmaz-hatsigmau)+#
 hatsigmazu=round(hatsigmaz-hatsigmau,8)
 chatsigma=(is.positive.definite(hatsigmazu))*hatsigmazu                 
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
 pvalue=rep(0,M)
 mm=0
 while(mm<M)
  {
  # Generating reproduce measurement error data tiaotaZ
   
   tiaotaZm=matrix(rep(0,n*d),n,d)     ### container of mean of tiaotaZ's at each hatX1
   tiaotaZsd=matrix(rep(0,n*d*d),n,d*d) ### container of covariance of tiaotaZ's at each hatX1
 
   for(i in 1:n)
    {
     Si=matrix(Zcov[i,],nrow=2)
     tiaotau=mvrnorm(nj,rep(0,d),Si)         ##produce measurement error tiaotau
     reptiaotaZ=t(matrix(rep(hatX1[i,],nj),2,nj))+tiaotau # repeat each row of hatX1 nj times
     tiaotaZm[i,]=apply(reptiaotaZ,2,mean)
 
     TiaotZ1=eigen(cov(reptiaotaZ))
     TiaotZ2=(TiaotZ1$vectors)%*%diag(sqrt(TiaotZ1$values))%*%t((TiaotZ1$vectors))
     tiaotaZsd[i,]=c(TiaotZ2)
    }
 
  # generating sample from the reproduce gamma mode regression 

   shapet=1+bett[d+3]
   ivegt=1/g(bett[1]+hatX1%*%bett[2:(d+1)]+X2*bett[d+2])
   scalet=bett[d+3]*ivegt
   tiaotay=rgamma(n,shape=shapet,rate=scalet)

   ZmBt=tiaotaZm[rep(seq_len(nrow(tiaotaZm)),each=B),] # repeat each row of tiaotaZm B times  
   ZvBt=tiaotaZsd[rep(seq_len(nrow(tiaotaZsd)),each =B), ] # repeat each row of tiaotaZcov B times 

  ZmBt1=ZmBt[1:(n*B/2),]
  ZmBt2=ZmBt[(n*B/2+1):(n*B),]

  ZvBt1=ZvBt[1:(n*B/2),]
  ZvBt2=ZvBt[(n*B/2+1):(n*B),]

  impartt=Impart(n,ZvBt)
  impartt1=Impart(n/2,ZvBt1)
  impartt2=Impart(n/2,ZvBt2)

  y1=tiaotay[1:(n/2)]
  y2=tiaotay[(n/2+1):n]
  

  mydata=list(y=tiaotay,data=data.frame(ZmB=ZmBt,impart=impartt,X2rep=X2rep))
  
  ests=optim(par=c(-0.3,-0.2,-0.2,-0.3,4),fn=mccl,mydata=mydata)
  

  
  if(ests$convergence==0)
   {
    point_estt = ests$par
    mm=mm+1
   } else 
   {
    next
   }
  

  mdata1=list(y=y1,data=data.frame(ZmB=ZmBt1,impart=impartt1,X2rep=X2rep1))
  mdata2=list(y=y2,data=data.frame(ZmB=ZmBt2,impart=impartt2,X2rep=X2rep2))

  
  Tvm=RTn(par=point_estt,mydata1=mdata1,mydata2=mdata2)
    pvalue[mm]=1*(Tvm>Tv)
    
    cat(mm, pvalue[mm], Tv, Tvm, "\n")
  }
  
  mean(pvalue)
  Pvalue=rbind(Pvalue,mean(pvalue))
  print(Pvalue)
  
  cat(jj,"\n")
}


mean(round(Pvalue,2)<0.05)
mean(round(Pvalue,2)<=0.05)

pv=seq(0,0.5,length=100)

fp1=fp2=rep(0,100)
for(j in seq(100))
{
fp1[j]=mean(round(Pvalue,2)<=pv[j])
fp2[j]=mean(round(Pvalue,2)<pv[j])
}

plot(pv,fp1,type="l",main="n=200.DC",xlab="Significance Level",ylab="Rejection Rate")
lines(pv,fp2)
#lines(pv,pv)


 te=Sys.time()
time=te-ts
save(file="n=200.M5dc.RData",list=ls())
