# Set working directory
setwd("//ad.ucl.ac.uk/homeu/rmjlmqu/Documents/Non-inferiority trials/NI Frontiers/Figures")

pdf(file="Base Case RD - Modify alpha.pdf",width=10,height=7)
m <- matrix(c(1,2,3,3),nrow = 2,ncol = 2,byrow = TRUE)
layout(mat = m,heights = c(0.75,0.25))

load("//ad.ucl.ac.uk/homeu/rmjlmqu/Documents/Non-inferiority trials/NI Frontiers/Simulation results/Base Case RD Modify alpha - T1Er.RData")
plot(range.of.pi0, T1Er.method[,1], type = "l", lwd=2, ylim=c(0,0.2),
     xlab = expression(pi[0]), ylab= "Type 1 Error Rate", main = "Risk difference",
     cex.lab=1.5, cex.axis=1.5, cex.main=2, cex.sub=1.5)
lines(range.of.pi0, T1Er.method[,2],  col="orange", type = "l", lwd=2, lty=2)
lines(range.of.pi0, T1Er.method[,3], col="red", type = "l", lwd=2, lty=4)
abline(h=0.025, col="red")

load("//ad.ucl.ac.uk/homeu/rmjlmqu/Documents/Non-inferiority trials/NI Frontiers/Simulation results/Base Case RD Modify alpha - Power.RData")
plot(range.of.pi0, Power.method[,1], type = "l", lwd=2, ylim=c(0,1), 
     xlab = expression(pi[0]), ylab= "Power", main = "Risk difference",
     cex.lab=1.5, cex.axis=1.5, cex.main=2, cex.sub=1.5)
lines(range.of.pi0, Power.method[,2],  col="orange", type = "l", lwd=2, lty=2)
lines(range.of.pi0, Power.method[,3], col="red", type = "l", lwd=2, lty=4)
abline(h=0.9, col="red")

par(mar = c(0.5,0.5,0.5,0.5))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
plot_colors <- c("black","orange", "red")
legend(x = "top",inset = 0,
       legend = c("Don't modify margin",expression("Modify margin - " ~ alpha ~ "=2.5%") ,expression("Modify margin - Choose" ~ alpha ~ "using Fig. 4 after observing " ~hat(pi)[0])), 
       col=plot_colors,  lwd=2, lty=c(1,2,4), cex=2)
par(mar = c(5.1, 4.1, 4.1, 2.1))
dev.off()