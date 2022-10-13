##################
#### Examples ####
##################

library(dani)

# Design parameters

n0<-n1<-166   # Number of observed outcomes in both arms
NI.marg<-0.1  # Original NI margin
sig.level.design<-0.025  # Original significance level
sig.level.modified<-0.005    # Lower significance level to maintain type 1 error under control when changing margin

# The SAFE frontier:
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


#Examples (comment out all but the one you want to do):

# events<-c(20, 20)  # Example 1: p0=p1=12%
events<-c(14,19)  # Example 2: p0=8.4%, p1=11.4%
# events<-c(2,9)  # Example 3: p0=1.2%, p1=5.4%

# call e0 and e1 event numbers in control and active and p0.unconstr the prob of event in control
e0<-events[1]
e1<-events[2]
p0.unconstr<-e0/n0

#Keep margin fixed:
NI.marg.kmf<-NI.marg
sig.level.kmf<-sig.level.design
test.kmf<-test.NI(n0,n1,e0,e1,NI.marg.kmf,sig.level.kmf)
CI.kmf<-test.kmf$CI
pvalue.kmf<-test.kmf$p


#Simply change margin:
NI.marg.scm<-NI.frontier(p0.unconstr)
sig.level.scm<-sig.level.design
test.scm<-test.NI(n0,n1,e0,e1,NI.marg.scm,sig.level.scm)
CI.scm<-test.scm$CI
pvalue.scm<-test.scm$p

#Always change significance level:
NI.marg.acsl<-NI.frontier(p0.unconstr)
sig.level.acsl<-sig.level.modified
test.acsl<-test.NI(n0,n1,e0,e1,NI.marg.acsl,sig.level.acsl)
CI.acsl<-test.acsl$CI
pvalue.acsl<-test.acsl$p

#Change sig level if change NI margin:
NI.marg.cslicm<-NI.frontier(e0/n0)
sig.level.cslicm<-ifelse(NI.frontier(p0.unconstr)==NI.marg,sig.level.design,
                      sig.level.modified)
test.cslicm<-test.NI(n0,n1,e0,e1,NI.marg.cslicm,sig.level.cslicm)
CI.cslicm<-test.cslicm$CI
pvalue.cslicm<-test.cslicm$p

# LRT:

p0.unconstr<-e0/n0
p1.unconstr<-e1/n1
if (p1.unconstr>p0.unconstr+NI.frontier(p0.unconstr)) p1.unconstr<-p0.unconstr+NI.frontier(p0.unconstr)

# COnstrained likelihood:
Lik.c<-function(p0) {
  l<-p0^e0*(1-p0)^(n0-e0)*(p0+NI.frontier(p0))^e1*(1-p0-NI.frontier(p0))^(n1-e1)
}

p0.constr<-optimize(Lik.c, c(0,1), maximum = TRUE)$maximum
p1.constr<-p0.constr+NI.frontier(p0.constr)

# Uncostrained likelihood
Lik.u<-function(p0,p1) {
  l<-p0^e0*(1-p0)^(n0-e0)*(p1)^e1*(1-p1)^(n1-e1)
}

max.lik.c<-Lik.u(p0.constr,p1.constr)
max.lik.u<-Lik.u(p0.unconstr,p1.unconstr)
ml.stat<-(-2*log(max.lik.c/max.lik.u))
pvalue.LRT<-(1-pchisq(ml.stat,1))

#Compute test-based CIs:
Z.value<-qnorm(pvalue.LRT)
sigma.corresponding<-abs((p1.unconstr-p0.unconstr-NI.frontier(p0.unconstr))/Z.value)
CI.LRT<-c(p1.unconstr-p0.unconstr-qnorm(0.975)*sigma.corresponding,
          p1.unconstr-p0.unconstr+qnorm(0.975)*sigma.corresponding)
sig.level.LRT<-sig.level.design
NI.marg.LRT<-NA

NI.margins<-c(NI.marg.kmf, NI.marg.scm, NI.marg.acsl, NI.marg.cslicm, NI.marg.LRT)
sig.levels<-c(sig.level.kmf, sig.level.scm, sig.level.acsl, sig.level.cslicm, sig.level.LRT)
lower.bounds<-c(CI.kmf[1], CI.scm[1], CI.acsl[1], CI.cslicm[1], CI.LRT[1])
upper.bounds<-c(CI.kmf[2], CI.scm[2], CI.acsl[2], CI.cslicm[2], CI.LRT[2])
pvalues<-c(pvalue.kmf, pvalue.scm, pvalue.acsl, pvalue.cslicm, pvalue.LRT)

results<-data.frame(NI.margins, sig.levels, lower.bounds, upper.bounds, pvalues)
row.names(results)<-c("Keep NI margin fixed", "Simply change NI margin",
                      "Reduced significance level", "Reduce when change NI margin",
                      "Likelihood Ratio Test")
colnames(results)<-c("NI margin", "Sig level", "Lower bound CI", "Upper bound CI", "p value")

View(results)