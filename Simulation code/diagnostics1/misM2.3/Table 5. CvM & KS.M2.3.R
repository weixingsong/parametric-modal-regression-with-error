# Simulation stuy for the paper "Parametric Modal Regression with Error Contaminated Covariates"
# Testing hypothesis: 
# Y~Gamma given link function and linear part are correctely specified.
#
# Test methods: Cramer von-Mise, and Kolmogrov-Smirnov    

 rm(list=ls())
 packages = c("MASS", "smoothmest","gofgamma", "VGAM")
 
 package.check = lapply( packages, FUN = function(x) 
    {
     if (!require(x, character.only = TRUE)) 
       {
         install.packages(x, dependencies = TRUE)
         library(x, character.only = TRUE)
       }
    } )

# mode-link function
 
 g=function(x){exp(x)}

# mccl function
 
 
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
 
#  Model Setup
 
 set.seed(4444)
 total=100
 
 pvseq=seq(0.01,0.5,by=0.01)
 lenpv=length(pvseq)
 
 rej1=rej2=matrix(0,nrow=total,ncol=lenpv)
 
 n=300   # sample size
 B=100   # B-size
 d=2     # X1,X2 are contaminated with normal error
 nj=4    # replication at each (X1,X2)
 mu=rep(0,d)   # mean of measurement error U
 su2=1      # variance of measurement error U
 covx1=diag(d)  # covariance of (X1,X2)
 covu=diag(rep(su2,d))  # covariance of u
 bt0=-0.25
 bt1=c(-0.25,-0.25)       # true beta-values for X1, X2
 bt2=-0.25                # true beta-value for X3
                
 
 jj=0
 while(jj<total)
  {
   phi=10         # true phi-value
   X2<-rbinom (n,1,0.5)    # generate  n sample from X3
   X1=matrix(0,nrow=n,ncol=2)# generate n sample from (X1,X2)
   for(i in seq(n))
    {
     X1[i,]<-mvrnorm(1,rep((X2[i]==1)-(X2[i]==0),d),covx1)
    }
  
   # generating sample from the mixture gamma mode regression 
 
   thetax=bt0+X1%*%bt1+bt2*X2
   alph1=1+phi               
   bet1=phi/g(thetax)
   bet2=phi/(1.3+g(thetax))
   un=runif(n,0,1)
   y=(un<0.3)*rgamma(n,shape=alph1,rate=bet1)+
     (un>0.3)*rgamma(n,shape=alph1,rate=bet2)
   
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

   X2rep=kronecker(X2,rep(1,B),"*")       # repeat each row of X2 B times
 
   mydata=list(y=y,data=data.frame(ZmB=ZmB,impart=impart,X2rep=X2rep))
   temp=optim(par=c(-0.3,-0.2,-0.2,-0.3,4),fn=mccl,mydata=mydata)
   
   if(temp$convergence==0)
     {
      bet = temp$par
      jj=jj+1
     } else {
       next
     }

    phi=bet[d+3]
    rp=ZmB%*%bet[2:(d+1)]+as.vector(bet[d+2]*X2rep)+bet[1]
    ip=impart%*%bet[2:(d+1)]
    Wb=complex(real=rp,imaginary=ip)
   
    invg=Re(apply(matrix(1/g(Wb),nrow=B),2,mean))
    R=y*phi*invg
  
    kk=1
    for(pv in pvseq) 
     {
      mytest1=test.CM(R, boot = 100, alpha = pv)$Decision
      mytest2=test.KS(R, boot = 100, alpha = pv)$Decision
      rej1[jj,kk]=mytest1
      rej2[jj,kk]=mytest2
      kk=kk+1
     }  
 
    cat(jj,bet,"\n")
  }

 p1=apply(rej1,2,mean)
 p2=apply(rej2,2,mean)
 ymaxi=max(p1,p2,pvseq)
 ymini=min(p1,p2,pvseq)
 #plot(pvseq,p1,type="l",col="blue",lwd=2,xlim=c(0.01,0.5),ylim=c(ymini,ymaxi))
 plot(pvseq,p1,type="l",col="blue",lwd=2,xlim=c(0.01,0.5),ylim=c(ymini,ymaxi),main="n=300.KS",xlab="Significance Level",ylab="Rejection Rate")
 lines(pvseq,p2,type="l",col="red",lwd=2)
# lines(pvseq,pvseq,lty=2,lwd=2)
 
#setwd("C:/Users/hp/Desktop/misM4")
#save(file="n=100.misM4.RData",list=ls())

save(file="n=300.ks0.5M4.RData",list=ls())