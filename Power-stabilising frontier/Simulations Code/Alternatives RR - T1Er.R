# Set seed and load libraries
set.seed(1)
library(epitools)

# Sample size calculation function:

sam_size = function(alpha,beta,pC,pA,NIm, r) {
  n= (qnorm(1-alpha)+qnorm(1-beta))^2*(((1-pC)/(pC)+(1-pA)/(r*pA))/(log(pA/pC)-NIm)^2)
  return(s<-c(n1<-ceiling(n),n2<- ceiling (ceiling(n)*r)))
}

# Parameterization alternatives

base.case.par<-list(
  n.sim<-100000,
  pie0<-0.05, # Expected control event rate
  pie1<-pie0, # Same as expected active event rate
  pif1<-0.1,  # Maximum tolerable active event rate 
  NI.marg<-log(pif1/pie0),
  r<-1,       # Allocation ratio
  power<-0.9, # Power
  alph<-0.025 # Significance level
)
names(base.case.par)<-c("n.sim", "pie0", "pie1", "pif1", "NI.marg", "r", "power", "alph")

n.alt<-7 # Number of alternatives

for (a in 1:n.alt) {
  assign(paste("alt",a,".par", sep=""),base.case.par)
}

alt1.par$pie0<-alt1.par$pie1<-0.1          # First alternative, different expecte dcontrol event rate
alt1.par$pif1<-2*alt1.par$pie0
alt1.par$NI.marg<-log(alt1.par$pif1/alt1.par$pie0)
alt2.par$pie1<-alt2.par$pie0/2             # Second alternative, expected active event rate is a half of control
alt3.par$pif1<-0.075                       # Third alternative, happy to tolerate up to 7.5%
alt3.par$NI.marg<-log(alt3.par$pif1/alt3.par$pie0)
alt4.par$pif1<-0.15                        # Fourth alternative, happy to tolerate up to 15%
alt4.par$NI.marg<-log(alt4.par$pif1/alt4.par$pie0)
alt5.par$r<-0.5                            # Fifth alternative, allocation ratio 0.5
alt6.par$r<-2                              # Sixth alternative, allocation ratio 2
alt7.par$power<-0.8                        # 80% power

# Initializations

range.of.pi0<-seq(0.005,0.2,0.005)
n.pi0<-length(range.of.pi0)
simulations<-array(NA, c(n.sim,42, n.pi0))
T1Er.method<-matrix(NA,n.pi0,7)
colnames(T1Er.method)<-1:7

# Iterate over different alternatives

for (a in 1:n.alt) {
  
  n.sim<-get(paste("alt",a,".par", sep=""))$n.sim
  pie0<-get(paste("alt",a,".par", sep=""))$pie0
  pie1<-get(paste("alt",a,".par", sep=""))$pie1
  pif1<-get(paste("alt",a,".par", sep=""))$pif1
  NI.marg<-get(paste("alt",a,".par", sep=""))$NI.marg
  r<-get(paste("alt",a,".par", sep=""))$r
  power<-get(paste("alt",a,".par", sep=""))$power
  alph<-get(paste("alt",a,".par", sep=""))$alph
  
  # Iterate over different values of "true" pi0
  
  for (j in 1:n.pi0) {
    pi0<-range.of.pi0[j]  # True control event rate
    pi1<-sin(asin(sqrt(pi0))+asin(sqrt(pif1))-asin(sqrt(pie0)))^2  # True active event rate
    # pi1<-ifelse(a==3,pi0*1.5,ifelse(a==4,pi0*3,pi0*2))  # True active event rate
    sample.size<-sam_size(alph, 1-power, pie0, pie1, NI.marg,r)  # Sample size
    
    # Data and results matrices
    
    simulations[,3+(a-1)*6,j]<-sample.size[1]
    simulations[,4+(a-1)*6,j]<-sample.size[2]
    
    for (i in 1:n.sim) {
      simulations[i,1+(a-1)*6,j]<-rbinom(1,sample.size[1],pi0)
      simulations[i,2+(a-1)*6,j]<-rbinom(1,sample.size[2],pi1)
      #simulations[i,5+(a-1)*4,1]<-prop.test(simulations[i,1:2,1],simulations[i,3:4,1],alternative = "g", correct = FALSE)$conf.int[1]
      # simulations[i,5+(a-1)*6,j]<-riskratio(matrix(c(simulations[i,1:2+(a-1)*6,j], simulations[i,3:4+(a-1)*6,j]-simulations[i,1:2+(a-1)*6,j]),ncol=2), rev = "c", conf.level = 0.9)$measure[2,3]
      selogRR<-sqrt(1/simulations[i,1+(a-1)*6,j]-1/sample.size[1]+1/simulations[i,2+(a-1)*6,j]-1/sample.size[2])
      logRR<-log((simulations[i,2+(a-1)*6,j]/sample.size[2])/(simulations[i,1+(a-1)*6,j]/sample.size[1]))
      simulations[i,5+(a-1)*6,j]<-exp(logRR+qnorm(1-alph)*selogRR)
      pio0<-simulations[i,1+(a-1)*6,j]/simulations[i,3+(a-1)*6,j]
      NI.marg2<-ifelse(abs(log(pio0/pie0))>log(1.5),
       #NI.marg2<-ifelse(abs(log(pio0/pie0))>log(1.25),
      # NI.marg2<-ifelse(abs(log(pio0/pie0))>log(2),
                       log(sin(asin(sqrt(pio0))+asin(sqrt(pif1))-asin(sqrt(pie0)))^2/pio0),
                       NI.marg)
      simulations[i,6+(a-1)*6,j]<-(log(simulations[i,5+(a-1)*6,j])<NI.marg2)
      
    }
    
    T1Er.method[j,a]<-sum(simulations[,6+(a-1)*6,j], na.rm=T)/(n.sim-sum(is.na(simulations[,6+(a-1)*6,j])))
    cat("Scenario ", j, "completed\n")
  }
  
  
}

# Save Results

# save.image("~/My Documents/Non-inferiority trials/NI Frontiers/Simulation results/Alternatives RR - T1Er.RData")
# save.image("~/My Documents/Non-inferiority trials/NI Frontiers/Simulation results/Alternatives 3 RR - T1Er.RData")
# save.image("~/My Documents/Non-inferiority trials/NI Frontiers/Simulation results/Alternatives 4 RR - T1Er.RData")
# save.image("~/My Documents/Non-inferiority trials/NI Frontiers/Simulation results/Alternatives RR - T1Er DGM Ian.RData")

# plot

par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=FALSE)
plot(range.of.pi0, T1Er.method[,1], type = "b", pch=16, ylim=c(0,0.4), xlab = expression(pi[0]), ylab= "Type 1 Error Rate", main = "Risk ratio")
lines(range.of.pi0, T1Er.method[,2],  col="blue", type = "b", pch=16)
lines(range.of.pi0, T1Er.method[,3], col="orange", type = "b", pch=16)
lines(range.of.pi0, T1Er.method[,4],  col="green", type = "b", pch=16, lty=2)
lines(range.of.pi0, T1Er.method[,5],  col="red", type = "b", pch=16, lty=2)
lines(range.of.pi0, T1Er.method[,6],  col="violet", type = "b", pch=16, lty=2)
lines(range.of.pi0, T1Er.method[,7],  col="grey", type = "b", pch=16, lty=2)
abline(h=0.05, col="red")
par(xpd=TRUE)
legend("topright",inset=c(-0.5,0),  legend = c(expression(pi[e0] ~ "=10%"), expression(pi[e1] ~ "=" ~ pi[e0]/2), expression(pi[f1] ~ "=7.5%"), expression(pi[f1] ~ "=15%"), "r=0.5", "r=2", "power=80%"), pch=16, col = c("black","blue","orange", "green", "red", "violet", "grey")) # optional legend
par(xpd=FALSE)
