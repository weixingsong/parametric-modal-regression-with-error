# Simulation stuy for the paper "Parametric Modal Regression with Error Contaminated Covariates"
# Testing hypothesis: 
# Y~Gamma given link function and linear part are correctely specified.
#
# Test methods: Cramer von-Mise, and Kolmogrov-Smirnov    

 rm(list=ls())

set.seed(6666)
# Load packages
 library(MASS)
 library(smoothmest)
 library(pracma) #complex number
 library(readxl)  # 载入readxl包
 library("scatterplot3d")
 library("plot3Drgl")
 library(matrixcalc)

# mode-link function
 
 g=function(x){exp(x)}
 d1g=function(x){exp(x)}


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
   Zsd=matrix(rep(0,d*d*n),nrow=n,ncol=d*d) # container of covariance of Z's at each X1
   Zcov=matrix(rep(0,d*d*n),nrow=n,ncol=d*d) # container of covariance of Z's at each X1
   for(j in seq(n))
    {
     Zmean[j,]=apply(Zrep[(nj*(j-1)+1):(nj*j),],2,mean)
     Zcov[j,]=c(cov(Zrep[(nj*(j-1)+1):(nj*j),]))

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

 Nresult=optim(c(-0.3,-0.2,-0.2,4),Naivmccl)$par
 
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
   
 
 result=optim(c(-0.3,-0.2,-0.2,4),mccl)$par
 
 An=hessian(mccl,result)
 Bn=mcov(result)
 invAn=solve(An)
 Bcov=sqrt(diag(invAn%*%Bn%*%invAn/n))

result
Nresult
Bcov
NBcov


####################HatX###############

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



#########################################

x1=hatX1[,1]
x2=hatX1[,2]

grid.lines = 200
x1.pred <- seq(min(x1), max(x1), length.out = grid.lines)
x2.pred <- seq(min(x2), max(x2), length.out = grid.lines)

x12 <- expand.grid( x = x1.pred, y = x2.pred)


z1pred= g(result[1]+result[2]*x12[,1]+result[3]*x12[,2])
z2pred= g(Nresult[1]+Nresult[2]*x12[,1]+Nresult[3]*x12[,2])


z1.pred <- matrix(z1pred, 
                 nrow = grid.lines, ncol = grid.lines)
z2.pred <- matrix(z2pred, 
                 nrow = grid.lines, ncol = grid.lines)


zlow=min(y)
zup=max(y)


# fitted points for droplines to surface
#fitpoints <- predict(fit)
# scatter plot with regression plane

par(mfrow=c(1,1))

x11()

scatter3D(x1, x2, y, pch = 19, cex = 1, colvar = NULL, col = "black",
    theta = 15, phi = 20, ticktype = "detailed",
    xlab = "Intake1", ylab = "Intake2", zlab = "FFQ",  
    surf = list(x = x1.pred, y = x2.pred, z = z1.pred,  
    facets = NA,
    col="gray66"), main = "FFQ vs. Long Time Intake(MCCL)",zlim=c(zlow,zup))


scatter3D(x1, x2, y, pch = 19, cex = 1, colvar = NULL, col = "black",
    theta = 15, phi = 20, ticktype = "detailed",
    xlab = "Intake1", ylab = "Intake2", zlab = "FFQ",  
    surf = list(x = x1.pred, y = x2.pred, z = z2.pred,  
    facets = NA,
    col="gray66"
    ), main = "FFQ vs. Long Time Intake(Naive)",zlim=c(zlow,zup),add=TRUE)
plotrgl()
