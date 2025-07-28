# Simulation stuy for the paper "Parametric Modal Regression with Error Contaminated Covariates"

rm(list=ls())
ts=Sys.time()

# Load packages
 library(MASS)
 library(smoothmest)
 library(pracma) #complex number

# mode-link function
 
 g=function(x){exp(x)}
 d1g=function(x){exp(x)}

 
set.seed(1234)

 total=1000 
 n=200   # sample size
 B=100   # B-size
 d=2     # X1,X2 are contaminated with normal error
 nj=4    # replication at each (X1,X2)
 mx1=rep(0,d)  # mean of (X1,X2)
 mu=rep(0,d)   # mean of measurement error U
 su2=1      # variance of measurement error U
 covx1=diag(d)  # covariance of (X1,X2)
 covu=diag(rep(su2,d))  # covariance of u

########################Naive#############

Nresult=NBcov=matrix(0,nrow=total,ncol=5) 

 for(kk in seq(total))
 {
 X2<-rbinom (n,1,0.5)    # generate  n sample from X3
 X1=matrix(0,nrow=n,ncol=2)# generate n sample from (X1,X2)
 for(i in seq(n))
  {
   X1[i,]<-mvrnorm(1,rep((X2[i]==1)-(X2[i]==0),d),covx1)
  }
 bt0=-0.25                # intercept
 bt1=c(-0.25,-0.25)       # true beta-values for X1, X2
 bt2=-0.25                # true beta-value for X3
 phi=5                    # true phi-value
 
 # generating sample from the gamma mode regression 
 
 alph=1+phi               
 bt=phi/g(bt0+X1%*%bt1+X2*bt2)
 y=rgamma(n,shape=alph,rate=bt)
 
 # Generating normal measurement error
 
 U=mvrnorm(nj*n,mu,covu)  
 X1rep=X1[rep(seq_len(nrow(X1)),each =nj), ]  # repeat each row of X1 nj times
 Zrep=X1rep+U   # nj replicated observations at each row of X1

 Zmean=matrix(rep(0,n*d),nrow=n)  # container of mean of Z's at each X1
 Zcov=matrix(rep(0,d*d*n),nrow=n,ncol=d*d) # container of covariance of Z's at each X1
 
 for(j in seq(n))
   {
    Zmean[j,]=apply(Zrep[(nj*(j-1)+1):(nj*j),],2,mean)
    # Spectral decomposition of covariance
    Temp1=eigen(cov(Zrep[(nj*(j-1)+1):(nj*j),]))
    Temp2=(Temp1$vectors)%*%diag(sqrt(Temp1$values))%*%t((Temp1$vectors))
    Zcov[j,]=c(Temp2)
   }

 Naivmccl=function(bet)
   {
    Wb=Zmean%*%bet[2:(d+1)]+as.vector(bet[d+2]*X2)+bet[1]
    logg=log(g(Wb))
    invg=1/g(Wb)
    phi=bet[d+3]
    out1=n*(1+phi)*log(phi)-n*log(gamma(phi+1))+phi*sum(log(y))
    out2=-(1+phi)*sum(logg)-phi*sum(y*invg)
    MCCL=out1+out2
    return(-MCCL/n)
   }
  
 Naivmcov=function(bet)
   {
    Wb01=Zmean[,1]
    Wb02=Zmean[,2]
    Wb=Zmean%*%bet[2:(d+1)]+as.vector(bet[d+2]*X2)+bet[1]
    phi=bet[d+3]
    
    Psi0= -(1+phi)*d1g(Wb)/g(Wb)+phi*y*d1g(Wb)/(g(Wb))^2
    Psi11= -(1+phi)*Wb01*d1g(Wb)/g(Wb)+phi*y*Wb01*d1g(Wb)/(g(Wb))^2
    Psi12= -(1+phi)*Wb02*d1g(Wb)/g(Wb)+phi*y*Wb02*d1g(Wb)/(g(Wb))^2
    Psi2= -(1+phi)*X2*d1g(Wb)/g(Wb)+phi*y*X2*d1g(Wb)/(g(Wb))^2
    Psi3= 1+log(phi)-digamma(phi)+log(y)-log(g(Wb))-y/g(Wb)
   
    Psi=cbind(Psi0,Psi11,Psi12,Psi2,Psi3)
    return((t(Psi)%*%Psi)/n)
   }

 Nresult[kk,]=optim(c(-0.3,-0.2,-0.2,-0.3,4),Naivmccl)$par
 
 An=hessian(Naivmccl,Nresult[kk,])
 Bn=Naivmcov(Nresult[kk,])
 invAn=solve(An)
 NBcov[kk,]=sqrt(diag(invAn%*%Bn%*%invAn/n))
 
 cat(kk,"\n")
 }

# boxplots for beta and phi
 par(mfrow=c(2,2))
 boxplot.matrix(Nresult[,1:4],xlab="Boxplots of Beta")
 abline(h=-0.25)
 boxplot(Nresult[,5],xlab="Boxplot of phi")
 abline(h=5)
 
 boxplot.matrix(NBcov[,1:4],xlab="Boxplots of std of Beta")
 boxplot(NBcov[,5],xlab="Boxplot of std of phi")


#Naive Numerical summary of beta and phi

 NBtM=apply(Nresult,2,mean)
 NBtMD=apply(Nresult,2,median)
 NBtS=apply(Nresult,2,sd)
 NBtIQ=apply(Nresult,2,IQR)

 NBtnumsum=rbind(NBtM,NBtMD,NBtS,NBtIQ)
 dimnames(NBtnumsum)=list(c("Mean","Median",
                         "Stdev","IQR"),
                       c("beta0","beta1","beta2","beta3","phi"))

#Naive Numerical summary of beta and phi  "sd"

 NSdM=apply(NBcov,2,mean)
 NSdMD=apply(NBcov,2,median)
 NSdS=apply(NBcov,2,sd)
 NSdIQ=apply(NBcov,2,IQR,na.rm=T)
 
 NSdnumsum=rbind(NSdM,NSdMD,NSdS,NSdIQ)
 dimnames(NSdnumsum)=list(c("Mean","Median",
                         "Stdev","IQR"),
                       c("beta0","beta1","beta2","beta3","phi"))

#######################MCCL###########
result=Bcov=matrix(0,nrow=total,ncol=5) 
# Generating Data
 for(jj in seq(total))
 {
 X2<-rbinom (n,1,0.5)    # generate  n sample from X3
 X1=matrix(0,nrow=n,ncol=2)# generate n sample from (X1,X2)
 for(i in seq(n))
  {
   X1[i,]<-mvrnorm(1,rep((X2[i]==1)-(X2[i]==0),d),covx1)
  }
 bt0=-0.25                # intercept
 bt1=c(-0.25,-0.25)       # true beta-values for X1, X2
 bt2=-0.25                # true beta-value for X3
 phi=5                    # true phi-value
 
 # generating sample from the gamma mode regression 
 
 alph=1+phi               
 bt=phi/g(bt0+X1%*%bt1+X2*bt2)
 y=rgamma(n,shape=alph,rate=bt)
 
 # Generating normal measurement error
 
 U=mvrnorm(nj*n,mu,covu)  
 X1rep=X1[rep(seq_len(nrow(X1)),each =nj), ]  # repeat each row of X1 nj times
 Zrep=X1rep+U   # nj replicated observations at each row of X1

 Zmean=matrix(rep(0,n*d),nrow=n)  # container of mean of Z's at each X1
 Zcov=matrix(rep(0,d*d*n),nrow=n,ncol=d*d) # container of covariance of Z's at each X1
 
 for(j in seq(n))
   {
    Zmean[j,]=apply(Zrep[(nj*(j-1)+1):(nj*j),],2,mean)
    # Spectral decomposition of covariance
    Temp1=eigen(cov(Zrep[(nj*(j-1)+1):(nj*j),]))
    Temp2=(Temp1$vectors)%*%diag(sqrt(Temp1$values))%*%t((Temp1$vectors))
    Zcov[j,]=c(Temp2)
   }
 
 ZmB=Zmean[rep(seq_len(nrow(Zmean)),each=B),] # repeat each row of Zmean B times  
 ZvB=Zcov[rep(seq_len(nrow(Zcov)),each =B), ] # repeat each row of Zcov B times 
 
 Tpm=mvrnorm(n*B,rep(0,nj-1),diag(nj-1))      # Generating Tp-values 
 Tpf=function(vv){vv[1:d]/sqrt(sum(vv^2))}
 Tp=t(apply(Tpm,1,Tpf))        
   
 Tp1=apply(ZvB[,1:d]*Tp,1,sum)          # Compute S*Tp
 Tp2=apply(ZvB[,(d+1):(2*d)]*Tp,1,sum)
 
 impart=sqrt((nj-1)/nj)*cbind(Tp1,Tp2)  # compute the imaginary part 
 realpart=ZmB                           # extract the real part
 
 X2rep=kronecker(X2,rep(1,B),"*")       # repeat each row of X2 B times
 
 
 #  Monte-Carlo Corrected Log-Likelihood 
 
 mccl=function(bet)
  {
   rp=ZmB%*%bet[2:(d+1)]+as.vector(bet[d+2]*X2rep)+bet[1]
   ip=impart%*%bet[2:(d+1)]
   Wb=complex(real=rp,imaginary=ip)
   logg=Re(apply(matrix(log(g(Wb)),nrow=B),2,mean))
   invg=Re(apply(matrix(1/g(Wb),nrow=B),2,mean))
   phi=bet[d+3]
   out1=n*(1+phi)*log(phi)-n*log(gamma(phi+1))+phi*sum(log(y))
   out2=-(1+phi)*sum(logg)-phi*sum(y*invg)
   MCCL=out1+out2
   return(-MCCL/n)
 }
 
 
 mcov=function(bet)
  {
    Wb01=complex(real=realpart[,1],imaginary=impart[,1])
    Wb02=complex(real=realpart[,2],imaginary=impart[,2])
    rp=ZmB%*%bet[2:(d+1)]+as.vector(bet[d+2]*X2rep)+bet[1]
    ip=impart%*%bet[2:(d+1)]
    Wb=complex(real=rp,imaginary=ip)
    phi=bet[d+3]
    yrep=kronecker(y,rep(1,B),"*")
    
    P0seq= -(1+phi)*d1g(Wb)/g(Wb)+phi*yrep*d1g(Wb)/(g(Wb))^2
    P1seq1= -(1+phi)*Wb01*d1g(Wb)/g(Wb)+phi*yrep*Wb01*d1g(Wb)/(g(Wb))^2
    P1seq2= -(1+phi)*Wb02*d1g(Wb)/g(Wb)+phi*yrep*Wb02*d1g(Wb)/(g(Wb))^2
    P2seq= -(1+phi)*X2rep*d1g(Wb)/g(Wb)+phi*yrep*X2rep*d1g(Wb)/(g(Wb))^2
    P3seq= 1+log(phi)-digamma(phi)+log(yrep)-log(g(Wb))-yrep/g(Wb)
    
    Psi0=Re(apply(matrix(P0seq,nrow=B),2,mean))
    Psi11=Re(apply(matrix(P1seq1,nrow=B),2,mean))
    Psi12=Re(apply(matrix(P1seq2,nrow=B),2,mean))
    Psi2=Re(apply(matrix(P2seq,nrow=B),2,mean))
    Psi3=Re(apply(matrix(P3seq,nrow=B),2,mean))
    Psi=rbind(Psi0,Psi11,Psi12,Psi2,Psi3)
    return((Psi%*%t(Psi))/n)
  }
   
 
 result[jj,]=optim(c(-0.3,-0.2,-0.2,-0.3,4),mccl)$par
 
 An=hessian(mccl,result[jj,])
 Bn=mcov(result[jj,])
 invAn=solve(An)
 Bcov[jj,]=sqrt(diag(invAn%*%Bn%*%invAn/n))
 
 cat(jj,"\n")
 }
 
 # boxplots for beta and phi


 par(mfrow=c(2,2))
 boxplot.matrix(result[,1:4],xlab="Boxplots of Beta")
 abline(h=-0.25)
 boxplot(result[,5],xlab="Boxplot of phi")
 abline(h=5)
 
 boxplot.matrix(Bcov[,1:4],xlab="Boxplots of std of Beta")
 boxplot(Bcov[,5],xlab="Boxplot of std of phi")
 
 
 
 # Numerical summary of beta and phi
 
 BtM=apply(result,2,mean)
 BtMD=apply(result,2,median)
 BtS=apply(result,2,sd)
 BtIQ=apply(result,2,IQR)

 Btnumsum=rbind(BtM,BtMD,BtS,BtIQ)
 dimnames(Btnumsum)=list(c("Mean","Median",
                         "Stdev","IQR"),
                       c("beta0","beta1","beta2","beta3","phi"))

 # Numerical summary of beta and phi
 
 SdM=apply(Bcov,2,mean)
 SdMD=apply(Bcov,2,median)
 SdS=apply(Bcov,2,sd)
 SdIQ=apply(Bcov,2,IQR)
 
 Sdnumsum=rbind(SdM,SdMD,SdS,SdIQ)
 dimnames(Sdnumsum)=list(c("Mean","Median",
                         "Stdev","IQR"),
                       c("beta0","beta1","beta2","beta3","phi"))




Btnumsum
Sdnumsum
NBtnumsum
NSdnumsum

resultNM=cbind(result[,1],Nresult[,1],result[,2],Nresult[,2],result[,3],Nresult[,3],result[,4],Nresult[,4])
boxplot.matrix(resultNM,xlab="Boxplots of Beta",main="X1 and X2 are dependent")
abline(h=-0.25)

te=Sys.time()
time=te-ts
save(file="Table2 sd 相依伯努利.RData",list=ls())

