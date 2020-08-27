# Set working directory
setwd("//ad.ucl.ac.uk/homeu/rmjlmqu/Documents/Non-inferiority trials/NI Frontiers/Figures")

pdf(file="Base Case - Modification rates.pdf",width=10,height=7)
m <- matrix(c(1,2,3,3),nrow = 2,ncol = 2,byrow = TRUE)
layout(mat = m,heights = c(0.75,0.25))

load("//ad.ucl.ac.uk/homeu/rmjlmqu/Documents/Non-inferiority trials/NI Frontiers/Simulation results/Base Case RD - T1Er.RData")
plot(range.of.pi0, n.modified[,1], type = "b", col="blue", pch=16, ylim=c(0,1), 
     xlab = expression(pi[0]), ylab= "Proportion modified", main = "Risk difference",
     cex.lab=1.5, cex.axis=1.5, cex.main=2, cex.sub=1.5)
lines(range.of.pi0, n.modified[,2], col="orange", type = "b", pch=16, lty=3)
lines(range.of.pi0, n.modified[,3],  col="green", type = "b", pch=16, lty=2)

load("//ad.ucl.ac.uk/homeu/rmjlmqu/Documents/Non-inferiority trials/NI Frontiers/Simulation results/Base Case RR - T1Er.RData")
plot(range.of.pi0, n.modified[,1], type = "b", col="blue", pch=16, ylim=c(0,1),
     xlab = expression(pi[0]), ylab= "Proportion modified", main = "Risk ratio",
     cex.lab=1.5, cex.axis=1.5, cex.main=2, cex.sub=1.5)
lines(range.of.pi0, n.modified[,2], col="orange", type = "b", pch=16, lty=3)
lines(range.of.pi0, n.modified[,3],  col="green", type = "b", pch=16, lty=2)

par(mar = c(0.5,0.5,0.5,0.5))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
plot_colors <- c("green","blue","orange")
legend(x = "top",inset = 0,
       legend = c(expression("Modify margin:" ~ epsilon ~ "=5% (RD) or log(2) (log-RR)"),expression("Modify margin:" ~ epsilon ~ "=2.5% (RD) or log(1.5) (log-RR)"), expression("Modify margin:" ~ epsilon ~ "=1.25% (RD) or log(1.25) (log-RR)")), 
       col=plot_colors,  pch=16, lty=c(2,1,3), cex=2)
par(mar = c(5.1, 4.1, 4.1, 2.1))
dev.off()