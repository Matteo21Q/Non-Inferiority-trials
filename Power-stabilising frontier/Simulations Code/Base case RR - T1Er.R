# Set seed and load libraries
set.seed(1)
library(epitools)

# Sample size calculation function:

sam_size = function(alpha,beta,pC,pA,NIm, r) {
  n= (qnorm(1-alpha)+qnorm(1-beta))^2*(((1-pC)/(pC)+(1-pA)/(r*pA))/(log(pA/pC)-NIm)^2)
  return(s<-c(n1<-ceiling(n),n2<- ceiling (ceiling(n)*r)))
}

# Parameterization

n.sim<-100000
pie0<-0.05 # Expected control event rate
pie1<-pie0 # Same as expected active event rate
pif1<-0.1  # Maximum tolerable active event rate 
NI.marg<-log(pif1/pie0)
r<-1       # Allocation ratio
power<-0.9 # Power
alph<-0.025 # Significance level

# Iterate over different values of "true" pi0

range.of.pi0<-seq(0.005,0.2,0.005)
n.pi0<-length(range.of.pi0)
simulations<-array(NA, c(n.sim,12, n.pi0))
T1Er.method<-matrix(NA,n.pi0,4)
colnames(T1Er.method)<-1:4
n.modified<-matrix(NA,n.pi0,3)
colnames(n.modified)<-2:4

for (j in 1:n.pi0) {
  pi0<-range.of.pi0[j]  # True control event rate
  pi1<-sin(asin(sqrt(pi0))+asin(sqrt(pif1))-asin(sqrt(pie0)))^2  # True active event rate
  # pi1<-pi0*2  # True active event rate
  sample.size<-sam_size(alph, 1-power, pie0, pie1, NI.marg,r)  # Sample size
  
  # Data and results matrices
  
  simulations[,3,j]<-sample.size[1]
  simulations[,4,j]<-sample.size[2]
  
  for (i in 1:n.sim) {
    simulations[i,1,j]<-rbinom(1,sample.size[1],pi0)
    simulations[i,2,j]<-rbinom(1,sample.size[2],pi1)
    #simulations[i,5,1]<-prop.test(simulations[i,1:2,1],simulations[i,3:4,1],alternative = "g", correct = FALSE)$conf.int[1]
    #simulations[i,5,j]<-riskratio(matrix(c(simulations[i,1:2,j], simulations[i,3:4,j]-simulations[i,1:2,j]),ncol=2), rev = "c", conf.level = 0.9)$measure[2,3]
    selogRR<-sqrt(1/simulations[i,1,j]-1/sample.size[1]+1/simulations[i,2,j]-1/sample.size[2])
    logRR<-log((simulations[i,2,j]/sample.size[2])/(simulations[i,1,j]/sample.size[1]))
    simulations[i,5,j]<-exp(logRR+qnorm(1-alph)*selogRR)
    simulations[i,6,j]<-(log(simulations[i,5,j])<NI.marg)
    pio0<-simulations[i,1,j]/simulations[i,3,j]
    NI.marg2<-ifelse(abs(log(pio0/pie0))>log(1.5),
                     log(sin(asin(sqrt(pio0))+asin(sqrt(pif1))-asin(sqrt(pie0)))^2/pio0),
                     NI.marg)
    simulations[i,7,j]<-(log(simulations[i,5,j])<NI.marg2)
    NI.marg3<-ifelse(abs(log(pio0/pie0))>log(1.25),
                     log(sin(asin(sqrt(pio0))+asin(sqrt(pif1))-asin(sqrt(pie0)))^2/pio0),
                     NI.marg)
    simulations[i,8,j]<-(log(simulations[i,5,j])<NI.marg3)
    NI.marg4<-ifelse(abs(log(pio0/pie0))>log(2),
                     log(sin(asin(sqrt(pio0))+asin(sqrt(pif1))-asin(sqrt(pie0)))^2/pio0),
                     NI.marg)
    simulations[i,9,j]<-(log(simulations[i,5,j])<NI.marg4)
    if (is.na(simulations[i,6:9,j])) simulations[i,6:9,j]<-0 # Too little events, hence equivalent to not being able to prove non-inferiority
    simulations[i,10,j]<-(NI.marg!=NI.marg2)
    simulations[i,11,j]<-(NI.marg!=NI.marg3)
    simulations[i,12,j]<-(NI.marg!=NI.marg4)
    
  }
  
  T1Er.method[j,1]<-sum(simulations[,6,j])/n.sim
  T1Er.method[j,2]<-sum(simulations[,7,j])/n.sim
  T1Er.method[j,3]<-sum(simulations[,8,j])/n.sim
  T1Er.method[j,4]<-sum(simulations[,9,j])/n.sim
  n.modified[j,1]<-sum(simulations[,10,j])/n.sim
  n.modified[j,2]<-sum(simulations[,11,j])/n.sim
  n.modified[j,3]<-sum(simulations[,12,j])/n.sim
  cat("Scenario ", j, "completed\n")
}

# save results

# save.image("~/My Documents/Non-inferiority trials/NI Frontiers/Simulation results/Base case RR - T1Er.RData")
# save.image("~/My Documents/Non-inferiority trials/NI Frontiers/Simulation results/Other DGM Results/Base case RR - T1Er DGM RR frontier.RData")

# Plot

par(mar=c(5.1, 4.1, 4.1, 6.1), xpd=FALSE)
plot(range.of.pi0, T1Er.method[,1], type = "b", pch=16, ylim=c(0,0.4), xlab = expression(pi[0]), ylab= "Type 1 Error Rate", main = "Riskratio")
lines(range.of.pi0, T1Er.method[,2],  col="blue", type = "b", pch=16)
lines(range.of.pi0, T1Er.method[,3], col="orange", type = "b", pch=16)
lines(range.of.pi0, T1Er.method[,4],  col="green", type = "b", pch=16, lty=2)
abline(h=0.05, col="red")
par(xpd=TRUE)
legend("topright",inset=c(-0.35,0), legend = c("Method 1", "Method 2", "Method 3", "Method 4"), pch=16, col = c("black","blue","orange", "green")) # optional legend
par(xpd=FALSE)

par(mar=c(5.1, 4.1, 4.1, 6.1), xpd=FALSE)
plot(range.of.pi0, n.modified[,1], type = "b",,  col="blue", pch=16, ylim=c(0,1), xlab = expression(pi[0]), ylab= "Proportion modified", main = "Risk ratio")
lines(range.of.pi0, n.modified[,2], col="orange", type = "b", pch=16)
lines(range.of.pi0, n.modified[,3],  col="green", type = "b", pch=16, lty=2)
par(xpd=TRUE)
legend("topright",inset=c(-0.35,0), legend = c("Method 2", "Method 3", "Method 4"), pch=16, col = c("blue","orange", "green")) # optional legend
par(xpd=FALSE)