
# Figure 7: Powers of MCCM1 in Model 2.1, 2.2, and 2.3

par(mfrow=c(1,3))

nsize=c(100,200,300)

model21=matrix(c(0.87,0.43,0.21,0.97,0.73,0.49,1.00,0.92,0.56),nrow=3,byrow=T)
plot(nsize,model21[,1],type="l",ylim=c(min(model21),max(model21)),
     xlab="n (Model 2.1)",ylab="Rejection Rate",lty=1,lwd=2,xaxt="n")
axis(side=1,at=c(100,200,300))
lines(nsize,model21[,2],lty=2,lwd=2)
lines(nsize,model21[,3],lty=3,lwd=2)
legend("bottomright",                     # position
       legend = c(expression(beta[4] == -0.25), expression(beta[4]==-0.08),expression(beta[4]==-0.05)), # legend text
       lwd = 2,                        # line width
       lty = c(1,2,3),
       y.intersp = 1.2,
       bty = "n")   

model22=matrix(c(0.97,0.38,0.29,0.26,1.00,0.70,0.69,0.52,1.00,0.78,0.82,0.54),nrow=3,byrow=T)
plot(nsize,model22[,1],type="l",ylim=c(min(model22),max(model22)),
     xlab="n (Model 2.2)",ylab="",lty=1,lwd=2,xaxt="n")
axis(side=1,at=c(100,200,300))
lines(nsize,model22[,2],lty=2,lwd=2)
lines(nsize,model22[,3],lty=3,lwd=2)
lines(nsize,model22[,4],lty=4,lwd=2)
legend("bottomright",                     # position
       legend = c(expression(beta[4] == 1), expression(beta[4]==0.27),expression(beta[4]==0.25),expression(beta[4]==0.2)), # legend text
       lwd = 2,                        # line width
       lty = c(1,2,3,4),
       y.intersp = 1.2,
       bty = "n")  

model23=matrix(c(0.73,0.57,0.47,0.40,0.98,0.85,0.76,0.58,1.00,0.90,0.84,0.70),nrow=3,byrow=T)
plot(nsize,model23[,1],type="l",ylim=c(min(model23),max(model23)),
     xlab="n (Model 2.3)",ylab="",lty=1,lwd=2,xaxt="n")
axis(side=1,at=c(100,200,300))
lines(nsize,model23[,2],lty=2,lwd=2)
lines(nsize,model23[,3],lty=3,lwd=2)
lines(nsize,model23[,4],lty=4,lwd=2)
legend("bottomright",                     # position
       legend = c(expression(beta[4] == 2), expression(beta[4]==1.4),expression(beta[4]==1.3),expression(beta[4]==1.2)), # legend text
       lwd = 2,                        # line width
       lty = c(1,2,3,4),
       y.intersp = 1.2,
       bty = "n") 

# Figure 8: Powers of MCCM1, MCCM2, DC, CvM and KS in Model 2.3

par(mfrow=c(1,1))
model=matrix(c(0.47,0.45,0.43,0.90,0.76,0.76,0.73,0.90,1.0,1.0,
               0.84,0.95,0.99,1,1),nrow=5,byrow=F)
plot(nsize,model[1,],type="l",ylim=c(min(model),max(model)),
     xlab="n (Model 2.3)",ylab="Rejection Rate",lty=1,lwd=2,xaxt="n")
axis(side=1,at=c(100,200,300))
lines(nsize,model[2,],lty=2,lwd=2)
lines(nsize,model[3,],lty=3,lwd=2)
lines(nsize,model[4,],lty=4,lwd=2)
lines(nsize,model[5,],lty=5,lwd=2)

legend("bottomright",                     # position
       legend = c("MCCM1","MCCM2","DC","CvM","KS"), # legend text
       lwd = 2,                        # line width
       lty = c(1,2,3,4,5),
       y.intersp = 1.2,
       bty = "n") 



