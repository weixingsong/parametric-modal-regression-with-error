# Simulation study for the paper "Parametric Modal Regression with Error Contaminated Covariates"
#  AD data analysis 
# Estimate parameters and standard (result：estimates of parameters，Bcov：standard errors).
#Note: mode-link function have three type as below.   The covariance of U have four type.


rm(list=ls())
ts=Sys.time()

# Load packages
 library(MASS)
 library(smoothmest)
 library(pracma) #complex number
 library(readxl)  # 载入readxl包

set.seed(6666)
 # mode-link function
 
 #g=function(x){exp(x)}
 #d1g=function(x){exp(x)}

 #g=function(x){1/(1+exp(x))}
 #d1g=function(x){-exp(x)/(1+exp(x))^2}

 g=function(x){exp(-exp(-x))}
 d1g=function(x){exp(-exp(-x))*exp(-x)}



 data<- read.csv("BJadni.csv", header=T)

 y=data$y
 X1=data[2:3]

 n=length(y)    # sample size
 B=50               # B-size
 d=2                 # X1 are contaminated with normal error
 nj=4                # replication at each (X1,X2)
 mu=rep(0,d)   # mean of measurement error U
 #covu=diag(c(0,0))          # covariance of u
 #covu=diag(c(0.03,0))     # covariance of u
 #covu=diag(c(0,0.16))     # covariance of u
 covu=diag(c(0.03,0.16))  # covariance of u


 U=mvrnorm(nj*n,mu,covu)  
 X1rep=X1[rep(seq_len(nrow(X1)),each =nj), ]  # repeat each row of X1 nj times
 Zrep=X1rep+U   # nj replicated observations at each row of X1

 Zmean=matrix(rep(0,n*d),nrow=n)                 # container of mean of Z's at each X1
 Zsd=matrix(rep(0,d*d*n),nrow=n,ncol=d*d)   # container of covariance of Z's at each X1
 
 for(j in seq(n))
   {
    Zmean[j,]=apply(Zrep[(nj*(j-1)+1):(nj*j),],2,mean)
    # Spectral decomposition of covariance
    Temp1=eigen(cov(Zrep[(nj*(j-1)+1):(nj*j),]))
    Temp2=(Temp1$vectors)%*%diag(sqrt(Temp1$values))%*%t((Temp1$vectors))
    Zsd[j,]=c(Temp2)
   }


########################Naive#############

 Naivmccl=function(bet)
   {
    Wb=Zmean%*%bet[2:(d+1)]+bet[1]
    logg=log(g(Wb))
    invg=1/g(Wb)
    phi=bet[d+2]
    out1=n*(1+phi)*log(phi)-n*log(gamma(phi+1))+phi*sum(log(y))
    out2=-(1+phi)*sum(logg)-phi*sum(y*invg)
    MCCL=out1+out2
    return(-MCCL/n)
   }

  

 Naivmcov=function(bet)
   {
    Wb01=Zmean[,1]
    Wb02=Zmean[,2]
    Wb=Zmean%*%bet[2:(d+1)]+bet[1]
    phi=bet[d+2]
    
    Psi0= -(1+phi)*d1g(Wb)/g(Wb)+phi*y*d1g(Wb)/(g(Wb))^2
    Psi11= -(1+phi)*Wb01*d1g(Wb)/g(Wb)+phi*y*Wb01*d1g(Wb)/(g(Wb))^2
    Psi12= -(1+phi)*Wb02*d1g(Wb)/g(Wb)+phi*y*Wb02*d1g(Wb)/(g(Wb))^2
    Psi3= 1+log(phi)-digamma(phi)+log(y)-log(g(Wb))-y/g(Wb)
   
    Psi=cbind(Psi0,Psi11,Psi12,Psi3)
    return((t(Psi)%*%Psi)/n)
   }

 Nresult=optim(c(-0.71,-0.15,-0.08,2.93),Naivmccl)$par
 
 An=hessian(Naivmccl,Nresult)
 Bn=Naivmcov(Nresult)
 invAn=solve(An)
 NBcov=sqrt(diag(invAn%*%Bn%*%invAn/n))
 

#######################MCCL###########

 ZmB=Zmean[rep(seq_len(nrow(Zmean)),each=B),] # repeat each row of Zmean B times  
 ZvB=Zsd[rep(seq_len(nrow(Zsd)),each =B), ] # repeat each row of Zcov B times 
 
 Tpm=mvrnorm(n*B,rep(0,nj-1),diag(nj-1))      # Generating Tp-values 
 Tpf=function(vv){vv[1:d]/sqrt(sum(vv^2))}
 Tp=t(apply(Tpm,1,Tpf))        
   
 Tp1=apply(ZvB[,1:d]*Tp,1,sum)          # Compute S*Tp
 Tp2=apply(ZvB[,(d+1):(2*d)]*Tp,1,sum)
 
 impart=sqrt((nj-1)/nj)*cbind(Tp1,Tp2)  # compute the imaginary part 
 realpart=ZmB                           # extract the real part
 
 
 
 #  Monte-Carlo Corrected Log-Likelihood 
 
 mccl=function(bet)
  {
   rp=ZmB%*%bet[2:(d+1)]+bet[1]
   ip=impart%*%bet[2:(d+1)]
   Wb=complex(real=rp,imaginary=ip)
   logg=Re(apply(matrix(log(g(Wb)),nrow=B),2,mean))
   invg=Re(apply(matrix(1/g(Wb),nrow=B),2,mean))
   phi=bet[d+2]
   out1=n*(1+phi)*log(phi)-n*log(gamma(phi+1))+phi*sum(log(y))
   out2=-(1+phi)*sum(logg)-phi*sum(y*invg)
   MCCL=out1+out2
   return(-MCCL/n)
  }
 

  
 
 mcov=function(bet)
  {
    Wb01=complex(real=realpart[,1],imaginary=impart[,1])
    Wb02=complex(real=realpart[,2],imaginary=impart[,2])
    rp=ZmB%*%bet[2:(d+1)]+bet[1]
    ip=impart%*%bet[2:(d+1)]
    Wb=complex(real=rp,imaginary=ip)
    phi=bet[d+2]
    yrep=kronecker(y,rep(1,B),"*")
    
    P0seq= -(1+phi)*d1g(Wb)/g(Wb)+phi*yrep*d1g(Wb)/(g(Wb))^2
    P1seq1= -(1+phi)*Wb01*d1g(Wb)/g(Wb)+phi*yrep*Wb01*d1g(Wb)/(g(Wb))^2
    P1seq2= -(1+phi)*Wb02*d1g(Wb)/g(Wb)+phi*yrep*Wb02*d1g(Wb)/(g(Wb))^2
    P3seq= 1+log(phi)-digamma(phi)+log(yrep)-log(g(Wb))-yrep/g(Wb)
    
    Psi0=Re(apply(matrix(P0seq,nrow=B),2,mean))
    Psi11=Re(apply(matrix(P1seq1,nrow=B),2,mean))
    Psi12=Re(apply(matrix(P1seq2,nrow=B),2,mean))
    Psi3=Re(apply(matrix(P3seq,nrow=B),2,mean))
    Psi=rbind(Psi0,Psi11,Psi12,Psi3)
    return((Psi%*%t(Psi))/n)
  }
   
 
 result=optim(c(-0.71,-0.15,-0.08,2.93),mccl)$par
 
 An=hessian(mccl,result)
 Bn=mcov(result)
 invAn=solve(An)
 Bcov=sqrt(diag(invAn%*%Bn%*%invAn/n))

result
Bcov



te=Sys.time()
time=te-ts

