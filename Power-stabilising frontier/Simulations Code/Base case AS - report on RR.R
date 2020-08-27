# Set seed and load libraries
set.seed(1)

# Sample size calculation function:

sam_size_RD= function(alpha,beta,pC,pA,NIm, r) {
  n= (qnorm(1-alpha)+qnorm(1-beta))^2*((pC*(1-pC)+pA*(1-pA)/r)/(pA-pC-NIm)^2)
  return(s<-c(n1<-ceiling(n),n2<- ceiling (ceiling(n)*r)))
}

sam_size_RR= function(alpha,beta,pC,pA,NIm, r) {
  n= (qnorm(1-alpha)+qnorm(1-beta))^2*(((1-pC)/(pC)+(1-pA)/(r*pA))/(log(pA/pC)-NIm)^2)
  return(s<-c(n1<-ceiling(n),n2<- ceiling (ceiling(n)*r)))
}

sam_size_AS= function(alpha,beta,pC,pA,NIm, r) {
  n= (qnorm(1-alpha)+qnorm(1-beta))^2*(((1/(4*r)+1/4))/(asin(sqrt(pA))-asin(sqrt(pC))-NIm)^2)
  return(s<-c(n1<-ceiling(n),n2<- ceiling (ceiling(n)*r)))
}

# Parameterization

n.sim<-10000
pie0<-0.05 # Expected control event rate
pie1<-pie0 # Same as expected active event rate
pif1<-0.1  # Maximum tolerable active event rate 
NI.marg<-asin(sqrt(pif1))-asin(sqrt(pie0))
NI.marg.RD<-pif1-pie0
NI.marg.RR<-log(pif1/pie0)
r<-1       # Allocation ratio
power<-0.9 # Power
alph<-0.025 # Significance level

# Iterate over different values of "true" pi0

range.of.pi0<-seq(0.005,0.2,0.005)
n.pi0<-length(range.of.pi0)
simulations<-array(NA, c(n.sim,8, n.pi0))
T1Er.method<-matrix(NA,n.pi0,4)
colnames(T1Er.method)<-1:4
n.modified<-matrix(NA,n.pi0,3)
colnames(n.modified)<-2:4
mean.margin<-rep(NA,n.pi0)
mean.alpha<-rep(NA,n.pi0)

for (j in 1:n.pi0) {
  pi0<-range.of.pi0[j]  # True control event rate
  pi1<-pi0  # True active event rate
  # pi1<-pi0*2  # True active event rate
  sample.size<-sam_size_AS(alph, 1-power, pie0, pie1, NI.marg,r)  # Sample size
  #sample.size<-sam_size_RD(alph, 1-power, pie0, pie1, NI.marg.RD,r)  # Sample size
  #sample.size<-sam_size_RR(alph, 1-power, pie0, pie1, NI.marg.RR,r)  # Sample size
  
  # Data and results matrices
  
  simulations[,3,j]<-sample.size[1]
  simulations[,4,j]<-sample.size[2]
  
  for (i in 1:n.sim) {
    simulations[i,1,j]<-rbinom(1,sample.size[1],pi0)
    simulations[i,2,j]<-rbinom(1,sample.size[2],pi1)
    #simulations[i,5,1]<-prop.test(simulations[i,1:2,1],simulations[i,3:4,1],alternative = "g", correct = FALSE)$conf.int[1]
    #simulations[i,5,j]<-riskratio(matrix(c(simulations[i,1:2,j], simulations[i,3:4,j]-simulations[i,1:2,j]),ncol=2), rev = "c", conf.level = 0.9)$measure[2,3]
    seAS<-sqrt(1/(4*sample.size[1])+1/(4*sample.size[2]))
    AS<-asin(sqrt(simulations[i,2,j]/sample.size[2]))-asin(sqrt(simulations[i,1,j]/sample.size[1]))
    simulations[i,5,j]<-AS+qnorm(1-alph)*seAS
    simulations[i,6,j]<-(simulations[i,5,j]<NI.marg)
    pio0<-simulations[i,1,j]/simulations[i,3,j]
    NI.marg.ASfront<-log(sin(asin(sqrt(pio0))+asin(sqrt(pif1))-asin(sqrt(pie0)))^2)-log(pio0)
    
    p.value<-pnorm((AS-NI.marg)/seAS)
    selogRR<-sqrt(1/simulations[i,1,j]-1/sample.size[1]+1/simulations[i,2,j]-1/sample.size[2])
    logRR<-log((simulations[i,2,j]/sample.size[2])/(simulations[i,1,j]/sample.size[1]))
    p.value.RR<-pnorm((logRR-NI.marg.RR)/selogRR)
    simulations[i,7,j]<-logRR-selogRR*qnorm(p.value)
    SE.RReq<-(logRR-NI.marg.ASfront)/qnorm(p.value)
    simulations[i,8,j]<-pnorm(qnorm(0.975)*SE.RReq/selogRR)
    
  }
  
  T1Er.method[j,1]<-sum(simulations[,6,j])/n.sim
  mean.margin[j]<-sum(simulations[,7,j])/n.sim
  mean.alpha[j]<-sum(simulations[,8,j])/n.sim
  
  cat("Scenario ", j, "completed\n")
}

yAS<-log(sin(asin(sqrt(range.of.pi0))+asin(sqrt(pif1))-asin(sqrt(pie0)))^2)-log(range.of.pi0)

# save results

# save.image("~/My Documents/Non-inferiority trials/NI Frontiers/Simulation results/Base case AS - Report on RR.RData")

# Plot


par(mar=c(5.1, 4.1, 4.1, 6.1), xpd=FALSE)
plot(range.of.pi0, mean.margin, type = "b", pch=16, ylim=c(0,3), xlab = expression(pi[0]), ylab= "NI margin", main = "Equivalent NI Margin")
lines(range.of.pi0, yAS,  col="blue", type = "b", pch=16)
abline(h=0.05, col="red")
par(xpd=TRUE)
legend("topright",inset=c(-0.35,0), legend = c("Method 1", "Method 2", "Method 3", "Method 4"), pch=16, col = c("black","blue","orange", "green")) # optional legend
par(xpd=FALSE)

par(mar=c(5.1, 4.1, 4.1, 6.1), xpd=FALSE)
plot(range.of.pi0, mean.alpha, type = "b", pch=16, ylim=c(0.9,1), xlab = expression(pi[0]), ylab= "1-sided level", main = "Equivalent confidence level")
abline(h=0.975, col="red")
par(xpd=TRUE)
legend("topright",inset=c(-0.35,0), legend = c("Method 1", "Method 2", "Method 3", "Method 4"), pch=16, col = c("black","blue","orange", "green")) # optional legend
par(xpd=FALSE)


