##########################
####  SAFE frontier  #####
####  Optimising alpha ###
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

# Modified significance level to try with the "modify sig level if modify margin" method
mod.sig.levs<-c(0.025,0.02,0.015,0.01,0.005)
# Initialise data.frame of results:
results.t1e<-data.frame(range.of.p0)

# Loop across significance levels:

for (sl in 1:length(mod.sig.levs)) {
  
  t1e.changeif<-NULL
  sig.level2<-mod.sig.levs[sl]
  
  for (p0 in range.of.p0) {
    p1 <- p0+NI.frontier(p0) # This is the only real diff: data generated at the null
    probs0<-rep(NA,sample.size[1]+1)
    probs1<-rep(NA,sample.size[2]+1)
    
    p.est.changeif<-0
    
    for (i in 0:sample.size[1]) {
      probs0[i+1]<-dbinom(i,sample.size[1],p0)
      if (probs0[i+1]>threshold.relevance) {
        pio0<-i/sample.size[1]
        
        NI.marg2 <- NI.frontier(pio0)
        
        p.est1.changeif<-0
        
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
          
          p.est1.changeif<-ifelse(pio0<0.09,p.est1.changeif+probs1[j+1]*(res.test.2<NI.marg2),
                                  p.est1.changeif+probs1[j+1]*(res.test.1<NI.marg1))
          
          
        }
        p.est.changeif<-p.est.changeif+probs0[i+1]*p.est1.changeif
        
      }
      
    }
    cat(".")
    t1e.changeif<-c(t1e.changeif,p.est.changeif)
    
  }
  results.t1e$newcol<-t1e.changeif
  colnames(results.t1e)[colnames(results.t1e)=="newcol"]<-paste("sig.lev.", mod.sig.levs[sl], sep = "")
  
  
}


##### View data frame of results:

View(results.t1e)
# Only 0.005 approximately controls type 1 error.
# Note that for low control risk t1e can be still ~3%, because we are using Wald intervals