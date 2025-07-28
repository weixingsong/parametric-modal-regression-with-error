# Simulation stuy for the paper  "Parametric Modal Regression with Error Contaminated Covariates"
# Testing hypothesis: 
# Y~Gamma given link function and linear part are correctely specified.
#
# Test methods: Cramer von-Mise, and Kolmogrov-Smirnov    

 rm(list=ls())
 packages = c("MASS", "smoothmest","pracma","readxl","gofgamma")
 
 package.check = lapply( packages, FUN = function(x) 
    {
     if (!require(x, character.only = TRUE)) 
       {
         install.packages(x, dependencies = TRUE)
         library(x, character.only = TRUE)
       }
    } )

set.seed(6666)
 
# mccl function
 
 
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



# mode-link function
 
 g=function(x){exp(x)}

 #g=function(x){1/(1+exp(x))}
 #d1g=function(x){-exp(x)/(1+exp(x))^2}

 #g=function(x){exp(-exp(-x))}
 # d1g=function(x){exp(-exp(-x))*exp(-x)}



 data<- read.csv("BJadni.csv", header=T)

 y=data$y
 X1=data[2:3]

 n=length(y)   # sample size
 B=50    # B-size
 d=2     # X1 are contaminated with normal error
 nj=4    # replication at each (X1,X2)
 mu=rep(0,d)   # mean of measurement error U
 #covu=diag(c(0,0))  # covariance of u
 #covu=diag(c(0.03,0))  # covariance of u
 #covu=diag(c(0,0.16))  # covariance of u
 covu=diag(c(0.03,0.16))  # covariance of u

 
   # Generating normal measurement error
 
   U=mvrnorm(nj*n,mu,covu)  
   X1rep=X1[rep(seq_len(nrow(X1)),each =nj), ]  # repeat each row of X1 nj times
   Zrep=X1rep+U   # nj replicated observations at each row of X1

   Zmean=matrix(rep(0,n*d),nrow=n)  # container of mean of Z's at each X1
   Zsd=matrix(rep(0,d*d*n),nrow=n,ncol=d*d) # container of covariance of Z's at each X1
 
   for(j in seq(n))
    {
     Zmean[j,]=apply(Zrep[(nj*(j-1)+1):(nj*j),],2,mean)
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
   bet=optim(par=c(-0.3,-0.2,-0.2,4),fn=mccl,mydata=mydata)$par

    phi=bet[d+2]
    rp=ZmB%*%bet[2:(d+1)]+bet[1]
    ip=impart%*%bet[2:(d+1)]
    Wb=complex(real=rp,imaginary=ip)
   
    invg=Re(apply(matrix(1/g(Wb),nrow=B),2,mean))
    R=y*phi*invg

      mytest1=test.CM(R, boot = 100, alpha = 0.05)$Decision
      mytest2=test.KS(R, boot = 100, alpha = 0.05)$Decision
      mytest1
      mytest2
      


