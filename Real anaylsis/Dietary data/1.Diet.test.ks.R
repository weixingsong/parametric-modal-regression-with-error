# Simulation stuy for the paper "Parametric Modal Regression with Error Contaminated Covariates" 
 rm(list=ls())
 packages = c("MASS", "smoothmest","pracma","matrixcalc","gofgamma")
 
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


    bet=betts$par
    phi=bet[3]
    rp=ZmB*bet[2]+bet[1]
    ip=impart*bet[2]
    Wb=complex(real=rp,imaginary=ip)
   
    invg=Re(apply(matrix(1/g(Wb),nrow=B),2,mean))
    R=y*phi*invg
  
   mytest1=test.CM(as.vector(R), boot = 100, alpha = 0.02)$Decision
   mytest2=test.KS(as.vector(R), boot = 100, alpha = 0.02)$Decision
   mytest1
   mytest2