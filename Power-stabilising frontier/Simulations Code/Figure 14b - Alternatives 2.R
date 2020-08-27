# Set working directory
setwd("~/Documents/Non-inferiority trials/NI Frontiers/Figures")

pdf(file="Alternatives 2.pdf",width=10,height=12)
m <- matrix(c(1,2,3,4,5,5),nrow = 3,ncol = 2,byrow = TRUE)
layout(mat = m,heights = c(0.4, 0.4, 0.2))

load("~/Documents/Non-inferiority trials/NI Frontiers/Simulation results/Alternatives RD - T1Er.RData")
plot(range.of.pi0, T1Er.method[,1], type = "l", lwd=2, ylim=c(0,0.2), 
     xlab = expression(pi[0]), ylab= "Type 1 Error Rate", main = "Risk difference",
     cex.lab=1.5, cex.axis=1.5, cex.main=2, cex.sub=1.5)
lines(range.of.pi0, T1Er.method[,2],  col="blue", type = "l", lwd=2, lty=2)
lines(range.of.pi0, T1Er.method[,3], col="orange", type = "l", lwd=4, lty=3)
lines(range.of.pi0, T1Er.method[,4],  col="green", type = "l", lwd=2, lty=4)
lines(range.of.pi0, T1Er.method[,5],  col="red", type = "l", lwd=2, lty=5)
lines(range.of.pi0, T1Er.method[,6],  col="violet", type = "l", lwd=2, lty=6)
lines(range.of.pi0, T1Er.method[,7],  col="grey", type = "l", lwd=2, lty=7)
abline(h=0.025, col="red")

load("~/Documents/Non-inferiority trials/NI Frontiers/Simulation results/Alternatives RR - T1Er.RData")
plot(range.of.pi0, T1Er.method[,1], type = "l", lwd=2, ylim=c(0,0.2), 
     xlab = expression(pi[0]), ylab= "Type 1 Error Rate", main = "Risk ratio",
     cex.lab=1.5, cex.axis=1.5, cex.main=2, cex.sub=1.5)
lines(range.of.pi0, T1Er.method[,2],  col="blue", type = "l", lwd=2, lty=2)
lines(range.of.pi0, T1Er.method[,3], col="orange", type = "l", lwd=4, lty=3)
lines(range.of.pi0, T1Er.method[,4],  col="green", type = "l", lwd=2, lty=4)
lines(range.of.pi0, T1Er.method[,5],  col="red", type = "l", lwd=2, lty=5)
lines(range.of.pi0, T1Er.method[,6],  col="violet", type = "l", lwd=2, lty=6)
lines(range.of.pi0, T1Er.method[,7],  col="grey", type = "l", lwd=2, lty=7)
abline(h=0.025, col="red")

load("~/DOcuments/Non-inferiority trials/NI Frontiers/Simulation results/Alternatives RD - Power.RData")
plot(range.of.pi0, Power.method[,1], type = "l", lwd=2, ylim=c(0,1), 
     xlab = expression(pi[0]), ylab= "Power", main = "Risk difference",
     cex.lab=1.5, cex.axis=1.5, cex.main=2, cex.sub=1.5)
lines(range.of.pi0, Power.method[,2],  col="blue", type = "l", lwd=2, lty=2)
lines(range.of.pi0, Power.method[,3], col="orange", type = "l", lwd=4, lty=3)
lines(range.of.pi0, Power.method[,4],  col="green", type = "l", lwd=2, lty=4)
lines(range.of.pi0, Power.method[,5],  col="red", type = "l", lwd=2, lty=5)
lines(range.of.pi0, Power.method[,6],  col="violet", type = "l", lwd=2, lty=6)
lines(range.of.pi0, Power.method[,7],  col="grey", type = "l", lwd=2, lty=7)

load("~/DOcuments/Non-inferiority trials/NI Frontiers/Simulation results/Alternatives RR - Power.RData")
plot(range.of.pi0, Power.method[,1], type = "l", lwd=2, ylim=c(0,1), 
     xlab = expression(pi[0]), ylab= "Power", main = "Risk ratio",
     cex.lab=1.5, cex.axis=1.5, cex.main=2, cex.sub=1.5)
lines(range.of.pi0, Power.method[,2],  col="blue", type = "l", lwd=2, lty=2)
lines(range.of.pi0, Power.method[,3], col="orange", type = "l", lwd=4, lty=3)
lines(range.of.pi0, Power.method[,4],  col="green", type = "l", lwd=2, lty=4)
lines(range.of.pi0, Power.method[,5],  col="red", type = "l", lwd=2, lty=5)
lines(range.of.pi0, Power.method[,6],  col="violet", type = "l", lwd=2, lty=6)
lines(range.of.pi0, Power.method[,7],  col="grey", type = "l", lwd=2, lty=7)

par(mar = c(0.5,0.5,0.5,0.5))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
plot_colors <- c("black","blue","orange", "green", "red", "violet", "grey")
legend(x = "top",inset = 0,
       legend = paste("Alternative ", 1:7), 
       col=plot_colors, lwd=c(2,2,4,2,2,2,2), lty=1:7, cex=2)
par(mar = c(5.1, 4.1, 4.1, 2.1))
dev.off()