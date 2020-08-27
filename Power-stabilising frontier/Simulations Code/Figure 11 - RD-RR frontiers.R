# Plot sample size given expected control event rate:

setwd("~/Documents/Non-inferiority trials/NI Frontiers/Figures")

pdf(file="NI Frontier.pdf",width=7,height=7)

m <- matrix(c(1,2),nrow = 2,ncol = 1,byrow = TRUE)
layout(mat = m,heights = c(0.85,0.15))

x<-seq(0,0.3,length.out = 1000)
yRD<-x+0.05
yRR<-2*x
yEQ<-x
yAS<-sin(asin(sqrt(x))+asin(sqrt(0.1))-asin(sqrt(0.05)))^2
yCI<-ifelse(x<0.025,x+0.025,ifelse(x<0.1,x+0.05,ifelse(x<0.2,x+0.075,x+0.1)))

plot(x,yEQ,type = "l", lwd=2, main="Non-inferiority frontiers (RD vs. RR vs. AS)", xlab = "Control event risk", ylab = "Active event risk")
lines(x,yRD, col="darkgoldenrod", lwd=2, lty=5)
lines(x,yRR,col="blue", lwd=2, lty=2)
lines(x,yAS,col="darkgreen", lwd=2, lty=4)
lines(x,yCI,col="grey", lwd=1)

par(mar = c(0.4,0.4,0.4,0.4))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
plot_colors <- c("black","darkgoldenrod","blue", "darkgreen", "grey")
legend(x = "top",inset = 0,
       legend = c("Treatment Equality Line","Fixed Risk Difference Frontier", "Fixed Risk Ratio Frontier", "Power-Stabilising Frontier", "Stepped Frontier"), 
       col=plot_colors, lty=c(1,5,2,4,1), cex=1, ncol=2, lwd=c(2,2,2,2,1))
par(mar = c(5.1, 4.1, 4.1, 2.1))

dev.off()