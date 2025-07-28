# Simulation stuy for the paper "Parametric Modal Regression with Error Contaminated Covariates"

rm(list=ls())
ts=Sys.time()

# Load packages
 library(MASS)
 library(smoothmest)
 library(pracma) #complex number

# mode-link function
 
 g=function(x){1/(1+exp(x))}
 d1g=function(x){-exp(x)/(1+exp(x))^2}

 
set.seed(1234)

 total=1000 
 n=2000   # sample size
 B=100   # B-size
 d=2     # X1,X2 are contaminated with normal error
 nj=4    # replication at each (X1,X2)
 mx1=rep(0,d)  # mean of (X1,X2)
 mu=rep(0,d)   # mean of measurement error U
 su2=1      # variance of measurement error U
 covx1=diag(d)  # covariance of (X1,X2)
 covu=diag(rep(su2,d))  # covariance of u

########################Naive#############

Nresult=matrix(0,nrow=total,ncol=5) 

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

  
 Nresult[kk,]=optim(c(-0.3,-0.2,-0.2,-0.3,4),Naivmccl)$par

 cat(kk,"\n") 
 }

# boxplots for beta and phi
 par(mfrow=c(2,2))
 boxplot.matrix(Nresult[,1:4],xlab="Boxplots of Beta")
 abline(h=-0.25)
 boxplot(Nresult[,5],xlab="Boxplot of phi")
 abline(h=5)

#Naive Numerical summary of beta and phi

 NBtM=apply(Nresult,2,mean)
 NBtMD=apply(Nresult,2,median)
 NBtS=apply(Nresult,2,sd)
 NBtIQ=apply(Nresult,2,IQR)

 NBtnumsum=rbind(NBtM,NBtMD,NBtS,NBtIQ)
 dimnames(NBtnumsum)=list(c("Mean","Median",
                         "Stdev","IQR"),
                       c("beta0","beta1","beta2","beta3","phi"))


#######################MCCL###########
result=matrix(0,nrow=total,ncol=5) 
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
 

  
 result[jj,]=optim(c(-0.3,-0.2,-0.2,-0.3,4),mccl)$par
 
 cat(jj,"\n")
 }
 
 # boxplots for beta and phi


 boxplot.matrix(result[,1:4],xlab="Boxplots of Beta")
 abline(h=-0.25)
 boxplot(result[,5],xlab="Boxplot of phi")
 abline(h=5) 
 
 
 # Numerical summary of beta and phi
 
 BtM=apply(result,2,mean)
 BtMD=apply(result,2,median)
 BtS=apply(result,2,sd)
 BtIQ=apply(result,2,IQR)

 Btnumsum=rbind(BtM,BtMD,BtS,BtIQ)
 dimnames(Btnumsum)=list(c("Mean","Median",
                         "Stdev","IQR"),
                       c("beta0","beta1","beta2","beta3","phi"))



Btnumsum
NBtnumsum


library(latex2exp)
resultNM=cbind(result[,1],Nresult[,1],result[,2],Nresult[,2],result[,3],Nresult[,3],result[,4],Nresult[,4])
boxplot.matrix(resultNM,xlab="MCCL(dark)    Naive(light)",col=c("gray61","gray94","gray61","gray94","gray61","gray94","gray61","gray94"),  at=c(1,2,4,5,7,8,10,11),xaxt="n")
xtick = c(1.5,4.5, 7.5, 10.5)
names=c(TeX(r"( $\beta_{0}$ )"),TeX(r"( $\beta_{1}$ )"),TeX(r"( $\beta_{2}$ )"),TeX(r"( $\beta_{3}$ )"))
axis(side=1, at=xtick, labels = names)
abline(h=-0.25,lty=3)
title(main = expression(bold(X)[1]~'a'*'n'*'d'~bold(X)[2]~'a'*'r'*'e'~'d'*'e'*'p'*'e'*'n'*'d'*'e'*'n'*'t'))


#resultNM=cbind(result[,1],Nresult[,1],result[,2],Nresult[,2],result[,3],Nresult[,3],result[,4],Nresult[,4])
#boxplot.matrix(resultNM,xlab="Boxplots of Beta",main="X1 and X2 are dependent")
#abline(h=-0.25)

te=Sys.time()
time=te-ts
save(file="figure1.相依伯努利.RData",list=ls())

