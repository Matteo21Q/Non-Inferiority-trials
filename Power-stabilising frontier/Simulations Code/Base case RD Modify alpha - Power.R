# Set seed and load libraries
set.seed(1)
library(Epi)

# Sample size calculation function:

sam_size = function(alpha,beta,pC,pA,NIm, r) {
  n= (qnorm(1-alpha)+qnorm(1-beta))^2*((pC*(1-pC)+pA*(1-pA)/r)/(pA-pC-NIm)^2)
  return(s<-c(n1<-ceiling(n),n2<- ceiling (ceiling(n)*r)))
}

# Parameterization

n.sim<-100000
pie0<-0.05 # Expected control event rate
pie1<-pie0 # Same as expected active event rate
pif1<-0.1  # Maximum tolerable active event rate 
NI.marg<-pif1-pie0
r<-1       # Allocation ratio
power<-0.9 # Power
alph<-0.025 # Significance level

# Iterate over different values of "true" pi0

range.of.pi0<-seq(0.005,0.2,0.005)
n.pi0<-length(range.of.pi0)
simulations<-array(NA, c(n.sim,9, n.pi0))
Power.method<-matrix(NA,n.pi0,3)
colnames(Power.method)<-1:3



for (j in 1:n.pi0) {
  pi0<-range.of.pi0[j]  # True control event rate
  pi1<-pi0  # True active event rate
  # pi1<-pi0+0.05  # True active event rate
  sample.size<-sam_size(alph, 1-power, pie0, pie1, -NI.marg,r)  # Sample size
  
  # Data and results matrices
  
  simulations[,3,j]<-sample.size[1]
  simulations[,4,j]<-sample.size[2]
  
  for (i in 1:n.sim) {
    simulations[i,1,j]<-rbinom(1,sample.size[1],pi0)
    simulations[i,2,j]<-rbinom(1,sample.size[2],pi1)
    #simulations[i,5,1]<-prop.test(simulations[i,1:2,1],simulations[i,3:4,1],alternative = "g", correct = FALSE)$conf.int[1]
    simulations[i,5,j]<-ci.pd(matrix(c(simulations[i,1:2,j], simulations[i,3:4,j]-simulations[i,1:2,j]),ncol=2, byrow=TRUE), alpha = 2* alph, print=F)[6]
    simulations[i,6,j]<-(abs(simulations[i,5,j])<NI.marg)
    pio0<-simulations[i,1,j]/simulations[i,3,j]
    NI.marg2<-ifelse(abs(pie0-pio0)>0.0125,
                     sin(asin(sqrt(pio0))+asin(sqrt(pif1))-asin(sqrt(pie0)))^2-pio0,
                     NI.marg)
    simulations[i,7,j]<-(abs(simulations[i,5,j])<NI.marg2)
    al.lev<-ifelse((pio0>0.005&&pio0<0.035),0.01,0.015)
    simulations[i,8,j]<-ci.pd(matrix(c(simulations[i,1:2,j], simulations[i,3:4,j]-simulations[i,1:2,j]),ncol=2, byrow=TRUE), alpha = al.lev*2, print=F)[6]
    simulations[i,9,j]<-(abs(simulations[i,8,j])<NI.marg2)
  }
  
  Power.method[j,1]<-sum(simulations[,6,j])/n.sim
  Power.method[j,2]<-sum(simulations[,7,j])/n.sim
  Power.method[j,3]<-sum(simulations[,9,j])/n.sim
  
  cat("Scenario ", j, "completed\n")
}

# Save results

# save.image("~/My Documents/Non-inferiority trials/NI Frontiers/Simulation results/Base case RD Modify alpha - Power.RData")

#Plot

par(mar=c(5.1, 4.1, 4.1, 6.1), xpd=FALSE)
plot(range.of.pi0, Power.method[,1], type = "b", pch=16, ylim=c(0,1), xlab = expression(pi[0]), ylab= "Power", main = "Risk difference")
lines(range.of.pi0, Power.method[,2],  col="orange", type = "b", pch=16)
lines(range.of.pi0, Power.method[,3], col="red", type = "b", pch=16)
abline(h=0.9, col="red")
par(xpd=TRUE)
legend("bottomleft",inset=c(-0.35,0), legend = c("Don't modify", "Modify margin but not alpha", "Modify margin  and alpha"), pch=16, col = c("black","orange", "red", "purple", "brown")) # optional legend
par(xpd=FALSE)
