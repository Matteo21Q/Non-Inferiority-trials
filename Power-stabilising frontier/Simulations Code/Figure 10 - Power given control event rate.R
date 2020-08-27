# Plot sample size given expected control event risk:

setwd("~/Documents/Non-inferiority trials/NI Frontiers/Figures")

pdf(file="Power - Control event risk.pdf",width=10,height=5)
m <- matrix(c(1,2),nrow = 1,ncol = 2,byrow = TRUE)
layout(mat = m,heights = c(0.85,0.85))

# Sample size calculation function for RD:

sam_size= function(alpha,beta,pC,pA,NIm, r) {
  n= (qnorm(1-alpha)+qnorm(1-beta))^2*((r*pC*(1-pC)+pA*(1-pA))/(pC-pA-NIm)^2)
  return(s<-c(n1<-ceiling(n),n2<- ceiling (ceiling(n)/r)))
}

alpha<-0.025
power<-0.9
pC<-seq(0.01,0.5,length.out = 100)
pA<-pC
NIm<-0.05
r<-1
ss<-sam_size(alpha,1-power,0.05,0.05,NIm,r)[1]

# Power calculation for RD

pow= function(alpha,samsize,pC,pA,NIm, r) {
  beta= (sqrt((samsize/((r*pC*(1-pC)+pA*(1-pA))/(pA-pC-NIm)^2)))-qnorm(1-alpha))
  return(p<-pnorm(beta))
}

powers<-rep(NA,length(pC))

for (i in 1:length(pC)) {
  powers[i]<-pow(alpha,ss,pC[i],pA[i],NIm,r)
}

plot(pC,powers,type = "l", xlim=c(0,0.3), las=1, lwd=2, ylab = "Power", xlab = "Control event risk", main="Risk difference (NI margin = 5%)")
abline(v=0.05, col="red")
abline(h=0.9, col="red")

# Sample size calculation function for RR:

sam_size_RR= function(alpha,beta,pC,pA,NIm, r) {
  n= (qnorm(1-alpha)+qnorm(1-beta))^2*((1/(r*pC)+1/(pA))/(log(pA/pC)-NIm)^2)
  return(s<-c(n1<-ceiling(n),n2<- ceiling (ceiling(n)*r)))
}

NImRR<-log(2)
ss2<-sam_size_RR(alpha,1-power,0.05,0.05,NImRR,r)[1]

# Power calculation for RR

pow2= function(alpha,samsize,pC,pA,NImRR, r) {
  beta= (sqrt((samsize/((1/(r*pC)+1/(pA))/(log(pA/pC)-NImRR)^2)))-qnorm(1-alpha))
  return(p<-pnorm(beta))
}

powers2<-rep(NA,length(pC))

for (i in 1:length(pC)) {
  powers2[i]<-pow2(alpha,ss2,pC[i],pA[i],NImRR,r)
}

plot(pC,powers2,type = "l", xlim=c(0,0.3), las=1, lwd=2, ylab = "Power", xlab = "Control event risk", main="Risk ratio (NI margin = 2)")
abline(v=0.05, col="red")
abline(h=0.9, col="red")

dev.off()