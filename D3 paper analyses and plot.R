##########################
####  SAFE frontier  #####
####  02-02-2022     #####
##########################

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
NI.frontier<-function(x) {
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
 
#Initialise power for various methods:
power.fixed<-NULL
power.always<-NULL
power.reduced<-NULL
power.changeif<-NULL
power.lrt<-NULL
 
# We loop for all possible values of true control risk:
for (p0 in range.of.p0) {
    p1 <- p0       # For power, assume same prob in control and active
    probs0<-rep(NA,sample.size[1]+1)
    probs1<-rep(NA,sample.size[2]+1)
    p.est.fixed<-0
    p.est.always<-0
    p.est.reduced<-0
    p.est.changeif<-0
    p.est.lrt<-0
    
    #Now we loop over all possible values of number of events in control:
    for (i in 0:sample.size[1]) {
      probs0[i+1]<-dbinom(i,sample.size[1],p0)  # Proba to observe that many events if p0
      if (probs0[i+1]>threshold.relevance) {
        pio0<-i/sample.size[1]
        
        NI.marg2 <- NI.frontier(pio0) # NI margin at that p0
        sig.level2<-0.005             # reduced sig level
        
        p.est1.fixed<-0
        p.est1.always<-0
        p.est1.reduced<-0
        p.est1.changeif<-0
        p.est1.lrt<-0
        
        for (j in 0:sample.size[2]) {
          probs1[j+1]<-dbinom(j,sample.size[2],p1) # Proba to observe that many events if p1
          
          # First, calculate CI for 2 possible sig levels:
          res.test.1<-test.NI(n0 = sample.size[1], 
                            n1 = sample.size[2], e0 = i, e1 = j, NI.margin = NI.marg2, 
                            sig.level = sig.level1, 
                            print.out = FALSE)$CI[2]
          res.test.2<-test.NI(n0 = sample.size[1], 
                              n1 = sample.size[2], e0 = i, e1 = j, NI.margin = NI.marg2, 
                              sig.level = sig.level2, 
                              print.out = FALSE)$CI[2]
          # Now implement the LRT method:
          p0.unconstr<-i/sample.size[1]
          p1.unconstr<-j/sample.size[2]
          if (p1.unconstr>p0.unconstr+NI.frontier(p0.unconstr)) p1.unconstr<-p0.unconstr+NI.frontier(p0.unconstr)
          
          # COnstrained likelihood:
          Lik.c<-function(p0) {
            l<-p0^i*(1-p0)^(sample.size[1]-i)*(p0+NI.frontier(p0))^j*(1-p0-NI.frontier(p0))^(sample.size[2]-j)
          }
          
          p0.constr<-optimize(Lik.c, c(0,1), maximum = TRUE)$maximum
          p1.constr<-p0.constr+NI.frontier(p0.constr)
          
          # Uncostrained likelihood
          Lik.u<-function(p0,p1) {
            l<-p0^i*(1-p0)^(sample.size[1]-i)*(p1)^j*(1-p1)^(sample.size[2]-j)
          }
          
          max.lik.c<-Lik.u(p0.constr,p1.constr)
          max.lik.u<-Lik.u(p0.unconstr,p1.unconstr)
          ml.stat<-(-2*log(max.lik.c/max.lik.u))
          p.value<-(1-pchisq(ml.stat,1)) # This will have to be halved

          p.est1.fixed<-p.est1.fixed+probs1[j+1]*(res.test.1<NI.marg1)
          p.est1.always<-p.est1.always+probs1[j+1]*(res.test.1<NI.marg2)
          p.est1.reduced<-p.est1.reduced+probs1[j+1]*(res.test.2<NI.marg2)
          p.est1.changeif<-ifelse(pio0<0.09,p.est1.changeif+probs1[j+1]*(res.test.2<NI.marg2),
                                  p.est1.changeif+probs1[j+1]*(res.test.1<NI.marg1))
          p.est1.lrt<-p.est1.lrt+probs1[j+1]*(p.value/2<sig.level1)
          
          
        }
        p.est.fixed<-p.est.fixed+probs0[i+1]*p.est1.fixed
        p.est.always<-p.est.always+probs0[i+1]*p.est1.always
        p.est.reduced<-p.est.reduced+probs0[i+1]*p.est1.reduced
        p.est.changeif<-p.est.changeif+probs0[i+1]*p.est1.changeif
        p.est.lrt<-p.est.lrt+probs0[i+1]*p.est1.lrt
        
      }
      
    }
     cat(".")
    power.fixed<-c(power.fixed,p.est.fixed)
    power.always<-c(power.always,p.est.always)
    power.reduced<-c(power.reduced,p.est.reduced)
    power.changeif<-c(power.changeif,p.est.changeif)
    power.lrt<-c(power.lrt,p.est.lrt)
    
  }
  results<-data.frame(range.of.p0,power.fixed, power.always, 
                      power.reduced, power.changeif, power.lrt)

  # Now repeat everythign but for type 1 error:
  t1e.fixed<-NULL
  t1e.always<-NULL
  t1e.reduced<-NULL
  t1e.changeif<-NULL
  t1e.lrt<-NULL
  
  for (p0 in range.of.p0) {
    p1 <- p0+NI.frontier(p0) # This is the only real diff: data generated at the null
    probs0<-rep(NA,sample.size[1]+1)
    probs1<-rep(NA,sample.size[2]+1)
    p.est.fixed<-0
    p.est.always<-0
    p.est.reduced<-0
    p.est.changeif<-0
    p.est.lrt<-0
    
    for (i in 0:sample.size[1]) {
      probs0[i+1]<-dbinom(i,sample.size[1],p0)
      if (probs0[i+1]>threshold.relevance) {
        pio0<-i/sample.size[1]
        
        NI.marg2 <- NI.frontier(pio0)
        sig.level2<-0.005
        
        p.est1.fixed<-0
        p.est1.always<-0
        p.est1.reduced<-0
        p.est1.changeif<-0
        p.est1.lrt<-0
        
        for (j in 0:sample.size[2]) {
          probs1[j+1]<-dbinom(j,sample.size[2],p1)
          
          res.test.1<-test.NI(n0 = sample.size[1], 
                              n1 = sample.size[2], e0 = i, e1 = j, NI.margin = NI.marg2, 
                              sig.level = sig.level1, 
                              print.out = FALSE)$CI[2]
          res.test.2<-test.NI(n0 = sample.size[1], 
                              n1 = sample.size[2], e0 = i, e1 = j, NI.margin = NI.marg2, 
                              sig.level = sig.level2, 
                              print.out = FALSE)$CI[2]
          p0.unconstr<-i/sample.size[1]
          p1.unconstr<-j/sample.size[2]
          if (p1.unconstr>p0.unconstr+NI.frontier(p0.unconstr)) p1.unconstr<-p0.unconstr+NI.frontier(p0.unconstr)
          
          Lik.c<-function(p0) {
            l<-p0^i*(1-p0)^(sample.size[1]-i)*(p0+NI.frontier(p0))^j*(1-p0-NI.frontier(p0))^(sample.size[2]-j)
          }
          
          p0.constr<-optimize(Lik.c, c(0,1), maximum = TRUE)$maximum
          p1.constr<-p0.constr+NI.frontier(p0.constr)
          
          Lik.u<-function(p0,p1) {
            l<-p0^i*(1-p0)^(sample.size[1]-i)*(p1)^j*(1-p1)^(sample.size[2]-j)
          }
          
          max.lik.c<-Lik.u(p0.constr,p1.constr)
          max.lik.u<-Lik.u(p0.unconstr,p1.unconstr)
          ml.stat<-(-2*log(max.lik.c/max.lik.u))
          p.value<-(1-pchisq(ml.stat,1))
          
          p.est1.fixed<-p.est1.fixed+probs1[j+1]*(res.test.1<NI.marg1)
          p.est1.always<-p.est1.always+probs1[j+1]*(res.test.1<NI.marg2)
          p.est1.reduced<-p.est1.reduced+probs1[j+1]*(res.test.2<NI.marg2)
          p.est1.changeif<-ifelse(pio0<0.09,p.est1.changeif+probs1[j+1]*(res.test.2<NI.marg2),
                                  p.est1.changeif+probs1[j+1]*(res.test.1<NI.marg1))
          p.est1.lrt<-p.est1.lrt+probs1[j+1]*(p.value/2<sig.level1)
          
          
        }
        p.est.fixed<-p.est.fixed+probs0[i+1]*p.est1.fixed
        p.est.always<-p.est.always+probs0[i+1]*p.est1.always
        p.est.reduced<-p.est.reduced+probs0[i+1]*p.est1.reduced
        p.est.changeif<-p.est.changeif+probs0[i+1]*p.est1.changeif
        p.est.lrt<-p.est.lrt+probs0[i+1]*p.est1.lrt
        
      }
      
    }
    cat(".")
    t1e.fixed<-c(t1e.fixed,p.est.fixed)
    t1e.always<-c(t1e.always,p.est.always)
    t1e.reduced<-c(t1e.reduced,p.est.reduced)
    t1e.changeif<-c(t1e.changeif,p.est.changeif)
    t1e.lrt<-c(t1e.lrt,p.est.lrt)
    
  }
  results.t1e<-data.frame(range.of.p0,t1e.fixed, t1e.always, 
                      t1e.reduced, t1e.changeif, t1e.lrt)
  
  
  # Now lot power and type 1 error:
  
  png("Power and t1e - D3.png", height = 480, width = 960)
  m <- matrix(c(1,2,3,3),nrow = 2,ncol = 2,byrow = TRUE)
  layout(mat = m,heights = c(0.78,0.22))
  plot(range.of.p0, results[,2], type="l", main = "Power comparison",
       xlab = "True control event risk", ylab = "Power", ylim = c(0.5,1),
       lwd=1)
  lines(range.of.p0,results[,3], col="red", lty=2)
  lines(range.of.p0,results[,4], col="blue", lty=3)
  lines(range.of.p0,results[,5], col="darkgreen", lty=4)
  lines(range.of.p0,results[,6], col="violet", lty=5)
  abline(h=0.8, col="grey")
  
  plot(range.of.p0, results.t1e[,2], type="l", main = "Type 1 error rate comparison",
       xlab = "True control event risk", ylab = "Type 1 error rate", ylim = c(0,0.1),
       lwd=1)
  lines(range.of.p0,results.t1e[,3], col="red", lty=2)
  lines(range.of.p0,results.t1e[,4], col="blue", lty=3)
  lines(range.of.p0,results.t1e[,5], col="darkgreen", lty=4)
  lines(range.of.p0,results.t1e[,6], col="violet", lty=5)
  abline(h=0.025, col="grey")
  par(mar=c(0.5,0.5,0.5,0.5))
  plot(1, type = "n", axes=FALSE, xlab="", ylab="")
  legend(x = "top",inset = 0,
         legend = c("Keep NI margin fixed","Simply change NI margin",
                    "Reduced significance level","Change sig. level when change NI margin",
                    "LRT"), 
         col=c("black","red","blue","darkgreen","violet"), 
         lty=c(1,2,3,4,5), ncol=1)
  dev.off()

  ### Plot frontier  ####
  
  y.front<-rep(NA,length(range.of.p0))
  for ( i in 1:length(range.of.p0)) {
    y.front[i]<-NI.frontier(range.of.p0[i])
    
  }
  # 
  png(paste("SAFE Frontier - D3.png", sep=""))
  plot(range.of.p0, y.front, type="l", main = paste("SAFE Frontier"), 
       xlab = "Control event risk", ylab = "Non-inferiority margin", ylim = c(0,0.15),
       lwd=3)
  abline(h=NI.marg1, col="red")
  dev.off()
  
