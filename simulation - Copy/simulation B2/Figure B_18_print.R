# Figure B.18
# Left Panel

par(mfrow=c(1,2))

nsize=c(100,200,300)
model14=matrix(c(0.53,0.61,1,0.94,0.65,0.83,1,0.99,0.64,0.99,1,1),nrow=4,byrow=F)

plot(nsize,model14[1,],type="l",ylim=c(min(model14)-0.1,max(model14)),
     xlab="n",ylab="Rejection Rate",lty=1,lwd=2,xaxt="n")
axis(side=1,at=c(100,200,300))
lines(nsize,model14[2,],lty=2,lwd=2)
lines(nsize,model14[3,],lty=3,lwd=2)
lines(nsize,model14[4,],lty=4,lwd=2)
legend("bottomright",                     # position
       legend = c("Model B.1","Model B.2","Model B.3","Model B.4"), # legend text
       lwd = 2,                        # line width
       lty = c(1,2,3,4),
       y.intersp = 1.2,
       bty = "n")   

# Table B.18
# Right Panel

model16=matrix(c(0.94,0.95,0.67,0.50,0.49,
                 0.99,1.00,0.94,0.86,0.87,
                 1.00,1.00,1.00,0.98,0.96
                 ),nrow=5,byrow=F)

model16
plot(nsize,model16[1,],type="l",ylim=c(min(model16),max(model16)),
     xlab="n",ylab="Rejection Rate",lty=1,lwd=2,xaxt="n")
axis(side=1,at=c(100,200,300))
lines(nsize,model16[2,],lty=2,lwd=2)
lines(nsize,model16[3,],lty=3,lwd=2)
lines(nsize,model16[4,],lty=4,lwd=2)
lines(nsize,model16[5,],lty=5,lwd=2)
legend("bottomright",                     # position
       legend = c("MCCM1","MCCM2","CvM","KS","DC"), # legend text
       lwd = 2,                        # line width
       lty = c(1,2,3,4,5),
       y.intersp = 1.2,
       bty = "n")   
