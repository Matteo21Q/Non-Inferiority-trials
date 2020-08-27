# Set working directory
setwd("~/Non-inferiority trials/NI Frontiers/Figures")

pdf(file="DGM.pdf",width=10,height=7)
m <- matrix(c(1,2,3,3),nrow = 2,ncol = 2,byrow = TRUE)
layout(mat = m,heights = c(0.85,0.15))

load("~/Non-inferiority trials/NI Frontiers/Simulation results/Base Case RD - T1Er.RData")
plot(range.of.pi0, T1Er.method[,1], type = "b", pch=16, ylim=c(0,0.4), xlab = expression(pi[0]), ylab= "Type 1 Error Rate", main = "Risk difference")
abline(h=0.05, col="red")
load("~/Non-inferiority trials/NI Frontiers/Simulation results/Base Case RD - T1Er DGM Ian.RData")
lines(range.of.pi0, T1Er.method[,1],  col="blue", type = "b", pch=16)

load("~/Non-inferiority trials/NI Frontiers/Simulation results/Base Case RR - T1Er.RData")
plot(range.of.pi0, T1Er.method[,1], type = "b", pch=16, ylim=c(0,0.4), xlab = expression(pi[0]), ylab= "Type 1 Error Rate", main = "Risk ratio")
abline(h=0.05, col="red")
load("~/Non-inferiority trials/NI Frontiers/Simulation results/Base Case RR - T1Er DGM Ian.RData")
lines(range.of.pi0, T1Er.method[,1],  col="blue", type = "b", pch=16)

par(mar = c(0.5,0.5,0.5,0.5))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
plot_colors <- c("black","blue")
legend(x = "top",inset = 0,
       legend = c("DGM I used", "DGM arcsin"), 
       col=plot_colors,  pch=16, cex=1, horiz = TRUE)
par(mar = c(5.1, 4.1, 4.1, 2.1))
dev.off()