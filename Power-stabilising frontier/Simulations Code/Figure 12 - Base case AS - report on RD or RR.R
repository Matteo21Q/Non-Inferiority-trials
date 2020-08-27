# Set working directory
setwd("~/My Documents/Non-inferiority trials/NI Frontiers/Figures")

pdf(file="Base Case AS - report on RD.pdf",width=10,height=7)
m <- matrix(c(1,2,3,4),nrow = 2,ncol = 2,byrow = TRUE)
layout(mat = m,heights = c(0.50,0.50))

load("~/My Documents/Non-inferiority trials/NI Frontiers/Simulation results/Base case AS - Report on RD.RData")
plot(range.of.pi0, mean.margin, type = "b", pch=16, ylim=c(0,0.15), xlab = expression(pi[0]), ylab= "NI margin", main = "Risk Difference")
lines(range.of.pi0, yAS,  col="blue", type = "b", pch=16)
abline(h=0.05, col="red")

plot(range.of.pi0, mean.alpha, type = "b", pch=16, ylim=c(0.9,1), xlab = expression(pi[0]), ylab= "1-sided level", main = "Risk Difference")
abline(h=0.975, col="red")

load("~/My Documents/Non-inferiority trials/NI Frontiers/Simulation results/Base case AS - Report on RR.RData")
plot(range.of.pi0, mean.margin, type = "b", pch=16, ylim=c(0,3), xlab = expression(pi[0]), ylab= "NI margin", main = "log-Risk Ratio")
lines(range.of.pi0, yAS,  col="blue", type = "b", pch=16)
abline(h=log(2), col="red")

plot(range.of.pi0, mean.alpha, type = "b", pch=16, ylim=c(0.9,1), xlab = expression(pi[0]), ylab= "1-sided level", main = "log-Risk Ratio")
abline(h=0.975, col="red")

dev.off()