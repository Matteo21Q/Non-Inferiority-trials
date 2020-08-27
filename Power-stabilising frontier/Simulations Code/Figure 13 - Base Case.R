# Set working directory
setwd("//ad.ucl.ac.uk/homeu/rmjlmqu/Documents/Non-inferiority trials/NI Frontiers/Figures")

pdf(file="Base Case.pdf",width=10,height=12)
m <- matrix(c(1,2,3,4,5,5),nrow = 3,ncol = 2,byrow = TRUE)
layout(mat = m,heights = c(0.425,0.425,0.15))

load("//ad.ucl.ac.uk/homeu/rmjlmqu/Documents/Non-inferiority trials/NI Frontiers/Simulation results/Base Case RD - T1Er.RData")
plot(range.of.pi0, T1Er.method[,1], type = "l", lwd=2, ylim=c(0,0.2), 
     xlab = expression(pi[0]), ylab= "Type 1 Error Rate", main = "Risk difference",
     cex.lab=1.5, cex.axis=1.5, cex.main=2, cex.sub=1.5)
lines(range.of.pi0, T1Er.method[,2],  col="blue", type = "l", lwd=4, lty=3)
lines(range.of.pi0, T1Er.method[,3], col="orange", type = "l", lwd=2, lty=4)
lines(range.of.pi0, T1Er.method[,4],  col="green", type = "l", lwd=2, lty=2)
abline(h=0.025, col="red")

load("//ad.ucl.ac.uk/homeu/rmjlmqu/Documents/Non-inferiority trials/NI Frontiers/Simulation results/Base Case RR - T1Er.RData")
plot(range.of.pi0, T1Er.method[,1], type = "l", lwd=2, ylim=c(0,0.2),
     xlab = expression(pi[0]), ylab= "Type 1 Error Rate", main = "Risk ratio",
     cex.lab=1.5, cex.axis=1.5, cex.main=2, cex.sub=1.5)
lines(range.of.pi0, T1Er.method[,2],  col="blue", type = "l", lwd=4, lty=3)
lines(range.of.pi0, T1Er.method[,3], col="orange", type = "l", lwd=2, lty=4)
lines(range.of.pi0, T1Er.method[,4],  col="green", type = "l", lwd=2, lty=2)
abline(h=0.025, col="red")

load("//ad.ucl.ac.uk/homeu/rmjlmqu/Documents/Non-inferiority trials/NI Frontiers/Simulation results/Base Case RD - Power.RData")
plot(range.of.pi0, Power.method[,1], type = "l", lwd=2, ylim=c(0,1), 
     xlab = expression(pi[0]), ylab= "Power", main = "Risk difference",
     cex.lab=1.5, cex.axis=1.5, cex.main=2, cex.sub=1.5)
lines(range.of.pi0, Power.method[,2],  col="blue", type = "l", lwd=4, lty=3)
lines(range.of.pi0, Power.method[,3], col="orange", type = "l", lwd=2, lty=4)
lines(range.of.pi0, Power.method[,4],  col="green", type = "l", lwd=2, lty=2)
abline(h=0.9, col="red")

load("//ad.ucl.ac.uk/homeu/rmjlmqu/Documents/Non-inferiority trials/NI Frontiers/Simulation results/Base Case RR - Power.RData")
plot(range.of.pi0, Power.method[,1], type = "l", lwd=2, ylim=c(0,1),
     xlab = expression(pi[0]), ylab= "Power", main = "Risk ratio",
     cex.lab=1.5, cex.axis=1.5, cex.main=2, cex.sub=1.5)
lines(range.of.pi0, Power.method[,2],  col="blue", type = "l", lwd=4, lty=3)
lines(range.of.pi0, Power.method[,3], col="orange", type = "l", lwd=2, lty=4)
lines(range.of.pi0, Power.method[,4],  col="green", type = "l", lwd=2, lty=2)
abline(h=0.9, col="red")

par(mar = c(0.5,0.5,0.5,0.5))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
plot_colors <- c("black", "green","blue","orange")
legend(x = "top",inset = 0,
       legend = c("Don't modify margin", expression("Modify margin:" ~ epsilon ~ "=5% (RD) or log(2) (log-RR)"),expression("Modify margin:" ~ epsilon ~ "=2.5% (RD) or log(1.5) (log-RR)"), expression("Modify margin:" ~ epsilon ~ "=1.25% (RD) or log(1.25) (log-RR)")), 
       col=plot_colors, lwd=c(2,2,4,2), lty=c(1,2,3,4), ncol=1, cex=2)
par(mar = c(5.1, 4.1, 4.1, 2.1))
dev.off()