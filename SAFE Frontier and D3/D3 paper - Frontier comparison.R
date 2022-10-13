################################
####  Frontier comparison  #####
####  02-02-2022           #####
################################

# Load relevant library 
library(dani)   # GitHub version 0.3.1

# Initialise parameters

threshold.relevance=10^(-5) # To decide when probability is non-negligible
range.of.p0<-seq(0.01,0.15,0.01) # Range of probs to consider
p0.expected<-0.12               # Expected event risk in control arm
NI.marg1<-0.1                   # Original NI margin   
sig.level1<-0.025               # Original one-sided alpha
power<-0.8                      # Original power level
sample.size<-sample.size.NI(p0.expected=p0.expected, 
                            p1.expected=p0.expected,
                            NI.margin=NI.marg1, 
                            sig.level = sig.level1,
                            power = power) 
sample.size<-10*signif(sample.size/10,2) # Sample size without increase for attrition

# Function to implement SAFE frontier:
NI.frontier1<-function(x) {
  NI.marg<-0.1
  threshold.difference=c(0.03)
  if (x>=0.09) {
    y<-NI.marg
  } else if ((x<0.09)&(x>0.05)) {
    y<-9/4*(x)-25/2*(x^2)-1/800
  } else {
    y<-0.75*x+0.0425
  }
  return(y)
}

# Function to implement power.stabilising frontier:
NI.frontier2<-function(x) {
  NI.marg<-0.1
  p0.expected<-0.12
  AS<-as.numeric(arcsine.margin(p0.expected+NI.marg,p0.expected))
  y<-optim(NI.marg,function(z) abs(arcsine.margin(z,x)-AS) , 
           method="L-BFGS-B", lower=0, upper=1)$par-x
  return(y)
}

# Function to implement Stepped frontier:
NI.frontier3<-function(x) {
  NI.marg<-0.1
  if (x>=0.09) {
    y<-NI.marg
  } else if ((x<0.09)&(x>=0.05)) {
    y<-0.08
  } else {
    y<-0.05
  }
  return(y)
}

# Function to implement Steep frontier:
NI.frontier4<-function(x) {
  NI.marg<-0.1
  if (x>=0.09) {
    y<-NI.marg
  } else if ((x<0.09)&(x>=0.08)) {
    y<-2*(x)-0.08
  } else if ((x<0.08)&(x>=0.05)) {
    y<-0.08
  } else if ((x<0.05)&(x>=0.04)) {
    y<-3*(x)-0.07
  } else {
    y<-0.05
  }
  return(y)
}

# Function to implement Linear frontier:
NI.frontier5<-function(x) {
  NI.marg<-0.1
  if (x>=0.09) {
    y<-NI.marg
  } else if ((x<0.09)&(x>0.05)) {
    y<-0.5*(x)+0.055
  } else {
    y<-0.75*x+0.0425
  }
  return(y)
}


#Initialise power for various frontiers:
power.SAFE<-NULL
power.arcsine<-NULL
power.steep<-NULL
power.stepped<-NULL
power.linear<-NULL

# We loop for all possible values of true control risk:
for (p0 in range.of.p0) {
  p1 <- p0       # For power, assume same prob in control and active
  probs0<-rep(NA,sample.size[1]+1)
  probs1<-rep(NA,sample.size[2]+1)
  p.est.SAFE<-0
  p.est.steep<-0
  p.est.stepped<-0
  p.est.linear<-0
  p.est.arcsine<-0
  
  #Now we loop over all possible values of number of events in control:
  for (i in 0:sample.size[1]) {
    probs0[i+1]<-dbinom(i,sample.size[1],p0)  # Proba to observe that many events if p0
    if (probs0[i+1]>threshold.relevance) {
      pio0<-i/sample.size[1]
      
      NI.marg1 <- NI.frontier1(pio0) # NI margin at that p0
      NI.marg2 <- NI.frontier2(pio0) # NI margin at that p0
      NI.marg3 <- NI.frontier3(pio0) # NI margin at that p0
      NI.marg4 <- NI.frontier4(pio0) # NI margin at that p0
      NI.marg5 <- NI.frontier5(pio0) # NI margin at that p0
      
      sig.level2<-sig.level1
      if (pio0<0.09) sig.level2<-0.005             # reduced sig level
      
      p.est1.SAFE<-0
      p.est1.steep<-0
      p.est1.stepped<-0
      p.est1.linear<-0
      p.est1.arcsine<-0
      
      for (j in 0:sample.size[2]) {
        probs1[j+1]<-dbinom(j,sample.size[2],p1) # Proba to observe that many events if p1
        
        # First, calculate CI for all possible frontiers:
        res.test.1<-test.NI(n0 = sample.size[1], 
                            n1 = sample.size[2], e0 = i, e1 = j, NI.margin = NI.marg1, 
                            sig.level = sig.level2, 
                            print.out = FALSE)$CI[2]
        res.test.2<-test.NI(n0 = sample.size[1], 
                            n1 = sample.size[2], e0 = i, e1 = j, NI.margin = NI.marg2, 
                            sig.level = 0.05, 
                            print.out = FALSE)$CI[2]
        res.test.3<-test.NI(n0 = sample.size[1], 
                            n1 = sample.size[2], e0 = i, e1 = j, NI.margin = NI.marg3, 
                            sig.level = sig.level2, 
                            print.out = FALSE)$CI[2]
        res.test.4<-test.NI(n0 = sample.size[1], 
                            n1 = sample.size[2], e0 = i, e1 = j, NI.margin = NI.marg4, 
                            sig.level = sig.level2, 
                            print.out = FALSE)$CI[2]
        res.test.5<-test.NI(n0 = sample.size[1], 
                            n1 = sample.size[2], e0 = i, e1 = j, NI.margin = NI.marg5, 
                            sig.level = sig.level2, 
                            print.out = FALSE)$CI[2]
        

        p.est1.SAFE<-p.est1.SAFE+probs1[j+1]*(res.test.1<NI.marg1)
        p.est1.arcsine<-p.est1.arcsine+probs1[j+1]*(res.test.2<NI.marg2)
        p.est1.steep<-p.est1.steep+probs1[j+1]*(res.test.3<NI.marg3)
        p.est1.stepped<-p.est1.stepped+probs1[j+1]*(res.test.4<NI.marg4)
        p.est1.linear<-p.est1.linear+probs1[j+1]*(res.test.5<NI.marg5)
        
        
      }
      p.est.SAFE<-p.est.SAFE+probs0[i+1]*p.est1.SAFE
      p.est.arcsine<-p.est.arcsine+probs0[i+1]*p.est1.arcsine
      p.est.steep<-p.est.steep+probs0[i+1]*p.est1.steep
      p.est.stepped<-p.est.stepped+probs0[i+1]*p.est1.stepped
      p.est.linear<-p.est.linear+probs0[i+1]*p.est1.linear
      
    }
    
  }
  cat(".")
  power.SAFE<-c(power.SAFE,p.est.SAFE)
  power.arcsine<-c(power.arcsine,p.est.arcsine)
  power.steep<-c(power.steep,p.est.steep)
  power.stepped<-c(power.stepped,p.est.stepped)
  power.linear<-c(power.linear,p.est.linear)
  
}
results<-data.frame(range.of.p0,power.SAFE, power.arcsine, 
                    power.steep, power.stepped, power.linear)

# Now repeat everythign but for type 1 error:
t1e.SAFE<-NULL
t1e.arcsine<-NULL
t1e.steep<-NULL
t1e.stepped<-NULL
t1e.linear<-NULL

for (p0 in range.of.p0) {
  p1.SAFE <- p0+NI.frontier1(p0) # This is the only real diff: data generated at the null
  p1.arcsine <- p0+NI.frontier2(p0) # This is the only real diff: data generated at the null
  p1.steep <- p0+NI.frontier3(p0) # This is the only real diff: data generated at the null
  p1.stepped <- p0+NI.frontier4(p0) # This is the only real diff: data generated at the null
  p1.linear <- p0+NI.frontier5(p0) # This is the only real diff: data generated at the null
  
  probs0<-rep(NA,sample.size[1]+1)
  probs1.SAFE<-rep(NA,sample.size[2]+1)
  probs1.arcsine<-rep(NA,sample.size[2]+1)
  probs1.steep<-rep(NA,sample.size[2]+1)
  probs1.stepped<-rep(NA,sample.size[2]+1)
  probs1.linear<-rep(NA,sample.size[2]+1)
  
  p.est.SAFE<-0
  p.est.arcsine<-0
  p.est.steep<-0
  p.est.stepped<-0
  p.est.linear<-0
  
  for (i in 0:sample.size[1]) {
    probs0[i+1]<-dbinom(i,sample.size[1],p0)
    if (probs0[i+1]>threshold.relevance) {
      pio0<-i/sample.size[1]
      
      NI.marg1 <- NI.frontier1(pio0) # NI margin at that p0
      NI.marg2 <- NI.frontier2(pio0) # NI margin at that p0
      NI.marg3 <- NI.frontier3(pio0) # NI margin at that p0
      NI.marg4 <- NI.frontier4(pio0) # NI margin at that p0
      NI.marg5 <- NI.frontier5(pio0) # NI margin at that p0
      
      sig.level2<-sig.level1
      if (pio0<0.09) sig.level2<-0.005             # reduced sig level
      
      p.est1.SAFE<-0
      p.est1.arcsine<-0
      p.est1.steep<-0
      p.est1.stepped<-0
      p.est1.linear<-0
      
      for (j in 0:sample.size[2]) {
        probs1.SAFE[j+1]<-dbinom(j,sample.size[2],p1.SAFE)
        probs1.arcsine[j+1]<-dbinom(j,sample.size[2],p1.arcsine)
        probs1.steep[j+1]<-dbinom(j,sample.size[2],p1.steep)
        probs1.stepped[j+1]<-dbinom(j,sample.size[2],p1.stepped)
        probs1.linear[j+1]<-dbinom(j,sample.size[2],p1.linear)
        
        # First, calculate CI for all possible frontiers:
        res.test.1<-test.NI(n0 = sample.size[1], 
                            n1 = sample.size[2], e0 = i, e1 = j, NI.margin = NI.marg1, 
                            sig.level = sig.level2, 
                            print.out = FALSE)$CI[2]
        res.test.2<-test.NI(n0 = sample.size[1], 
                            n1 = sample.size[2], e0 = i, e1 = j, NI.margin = NI.marg2, 
                            sig.level = 0.05, 
                            print.out = FALSE)$CI[2]
        res.test.3<-test.NI(n0 = sample.size[1], 
                            n1 = sample.size[2], e0 = i, e1 = j, NI.margin = NI.marg3, 
                            sig.level = sig.level2, 
                            print.out = FALSE)$CI[2]
        res.test.4<-test.NI(n0 = sample.size[1], 
                            n1 = sample.size[2], e0 = i, e1 = j, NI.margin = NI.marg4, 
                            sig.level = sig.level2, 
                            print.out = FALSE)$CI[2]
        res.test.5<-test.NI(n0 = sample.size[1], 
                            n1 = sample.size[2], e0 = i, e1 = j, NI.margin = NI.marg5, 
                            sig.level = sig.level2, 
                            print.out = FALSE)$CI[2]

        p.est1.SAFE<-p.est1.SAFE+probs1.SAFE[j+1]*(res.test.1<NI.marg1)
        p.est1.arcsine<-p.est1.arcsine+probs1.arcsine[j+1]*(res.test.2<NI.marg2)
        p.est1.steep<-p.est1.steep+probs1.steep[j+1]*(res.test.3<NI.marg3)
        p.est1.stepped<-p.est1.stepped+probs1.stepped[j+1]*(res.test.4<NI.marg4)
        p.est1.linear<-p.est1.linear+probs1.linear[j+1]*(res.test.5<NI.marg5)
        
        
      }
      p.est.SAFE<-p.est.SAFE+probs0[i+1]*p.est1.SAFE
      p.est.arcsine<-p.est.arcsine+probs0[i+1]*p.est1.arcsine
      p.est.steep<-p.est.steep+probs0[i+1]*p.est1.steep
      p.est.stepped<-p.est.stepped+probs0[i+1]*p.est1.stepped
      p.est.linear<-p.est.linear+probs0[i+1]*p.est1.linear
      
    }
    
  }
  cat(".")
  t1e.SAFE<-c(t1e.SAFE,p.est.SAFE)
  t1e.arcsine<-c(t1e.arcsine,p.est.arcsine)
  t1e.steep<-c(t1e.steep,p.est.steep)
  t1e.stepped<-c(t1e.stepped,p.est.stepped)
  t1e.linear<-c(t1e.linear,p.est.linear)
  
}
results.t1e<-data.frame(range.of.p0,t1e.SAFE, t1e.arcsine, 
                        t1e.steep, t1e.stepped, t1e.linear)


# Now plot power and type 1 error:

png("Frontier comparison - D3.png", height = 480, width = 960)
m <- matrix(c(1,2,3,3),nrow = 2,ncol = 2,byrow = TRUE)
layout(mat = m,heights = c(0.78,0.22))
plot(range.of.p0, results[,2], type="l", main = "Power comparison",
     xlab = "True control event risk", ylab = "Power", ylim = c(0.5,1),
     lwd=3)
lines(range.of.p0,results[,4], col="blue", lty=3)
lines(range.of.p0,results[,5], col="darkgreen", lty=4)
lines(range.of.p0,results[,6], col="violet", lty=5)
abline(h=0.8, col="grey")

plot(range.of.p0, results.t1e[,2], type="l", main = "Type 1 error rate comparison",
     xlab = "True control event risk", ylab = "Type 1 error rate", ylim = c(0,0.1),
     lwd=3)
lines(range.of.p0,results.t1e[,4], col="blue", lty=3)
lines(range.of.p0,results.t1e[,5], col="darkgreen", lty=4)
lines(range.of.p0,results.t1e[,6], col="violet", lty=5)
abline(h=0.025, col="grey")
par(mar=c(0.5,0.5,0.5,0.5))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
legend(x = "top",inset = 0,
       legend = c("Stepped frontier",
                  "Steep frontier",
                  "Linear frontier",
                  "SAFE frontier"),  
       col=c("darkgreen","blue","violet","black"), 
       lty=c(4,3,5,1), lwd=c(1,1,1,3), ncol=1)
dev.off()

### Plot frontiers  ####

y.front1<-y.front2<-y.front4<-y.front5<-rep(NA,length(range.of.p0))
range.of.p0.2<-c(0.01,0.02,0.03,0.04,0.0499,0.05,0.06,0.07,0.08,0.08999,0.09,
                 0.1,0.11,0.12,0.13,0.14,0.15)
y.front3<-rep(NA,length(range.of.p0.2))
for ( i in 1:length(range.of.p0)) {
  y.front1[i]<-NI.frontier1(range.of.p0[i])
  y.front2[i]<-NI.frontier2(range.of.p0[i])
  y.front4[i]<-NI.frontier4(range.of.p0[i])
  y.front5[i]<-NI.frontier5(range.of.p0[i])
}
for ( i in 1:length(range.of.p0.2)) {
  y.front3[i]<-NI.frontier3(range.of.p0.2[i])
}

# 
png(paste("All Frontiers - D3.png", sep=""), height = 600 )
m <- matrix(c(1,2),nrow = 2,ncol = 1,byrow = TRUE)
layout(mat = m,heights = c(0.78,0.22))
plot(range.of.p0.2, y.front3, type="l", main = paste("D3 Frontiers"), 
     xlab = "Control event risk", ylab = "Non-inferiority margin", ylim = c(0,0.15),
      col="darkgreen", lty=4)
lines(range.of.p0,y.front4, col="blue", lty=3)
lines(range.of.p0,y.front5, col="violet", lty=5)
lines(range.of.p0,lwd=3,y.front1)
abline(h=NI.marg1, col="red")
par(mar=c(0.5,0.5,0.5,0.5))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
legend(x = "top",inset = 0,
       legend = c("Stepped frontier",
                  "Steep frontier",
                  "Linear frontier",
                  "SAFE frontier",
                  "Fixed risk difference"), 
       col=c("darkgreen","blue","violet","black", "red"), 
       lty=c(4,3,5,1,1), lwd=c(1,1,1,3,1), ncol=1)
dev.off()
