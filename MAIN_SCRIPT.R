# set working directory to repo's location

# libraries needed
library(gdata) 
library(signal)
library(lme4)
library(pander)
library(stargazer)

source('functions.R') # variety of functions 

# load prior analyses for DVs
load('loadedSignals.RData') # the time series of body movements, from xls (see bodyExtract.R)

#
# source('computeCorrelations.R') # to build the RData files below
#

load('wccData.RData') # windowed correlation values, with surrogates
load('ccfresData.RData') # cross correlation functions across triads, with surrogates
load('wccfullData.RData') # the time series of windowed cross correlatiojs (for triadic)

wccres$max.loc = wccres$max.loc - 100 # center the lag term (1,201 to -100,100)

#
# triadic synchrony observed vs. surrogate
#

# 10s window analysis
lmo = lmer(winmax~cond+(1|triad),data=wccres) # max struct not possible
pander(print_stats(lmo))
hist(resid(lmo)) # check approx normal

lmo = lmer(winmin~cond+(1+cond|triad),data=wccres) # minimum observed simultaneous correlation
pander(print_stats(lmo))
hist(resid(lmo)) # check approx normal

# one-sample t-test to show that minimum correlation is, of course, significant
t.test(aggregate(winmin~triad,data=wccres[wccres$cond=='obs',],FUN=mean)$winmin)

# since some cannot converge, let's double check with paired samples t-test for winmin, winmax
agg=aggregate(winmin~triad+cond,data=wccres,mean) # switch winmax/winmin to check each
agg=cbind(agg[1:35,],agg[36:70,])
mean(agg[1:35,1]==agg[1:35,4]) # ensure alignment
t.test(agg[,3]-agg[,6])

# cross correlation function
lmo = lmer(ccf~cond+(1+cond|triad),data=ccfres[ccfres$lag==0,]) # correlation higher at lag 0
pander(print_stats(lmo))
hist(resid(lmo)) # check approx normal

# let's double check with paired samples t-test for ccf
agg=aggregate(ccf~triad+cond,data=ccfres[ccfres$lag==0,],mean)
agg=cbind(agg[1:35,],agg[36:70,])
mean(agg[1:35,1]==agg[1:35,4]) # ensure alignment
t.test(agg[,3]-agg[,6])

# reliable for each pair, too; using lag of 0 comparing to surrogate (vrt) and observed (obs)
lmo = lmer(ccf~cond+(1+cond|triad),data=ccfres[ccfres$lag==0&ccfres$typ=='LC',])
pander(print_stats(lmo))
lmo = lmer(ccf~cond+(1+cond|triad),data=ccfres[ccfres$lag==0&ccfres$typ=='CR',])
pander(print_stats(lmo))
lmo = lmer(ccf~cond+(1+cond|triad),data=ccfres[ccfres$lag==0&ccfres$typ=='LR',])
pander(print_stats(lmo))

#
# triad synchrony: is the triad itself moving together?
#

# let's get the CCFs of the correlations: OBSERVED (obs)
# this looks more complicated than it is... simply: loop through triads
# and for each compute the cross correlation of their windowed correlations (over time)
# since we have 3 people, we do this three times (3 pairs)... then average
# *NB: wccCcf = cross-correlation function of the windowed cross-correlations (second-order)
triads = unique(wccfull$triad)
for (i in 1:35) {
  print(i)
  
  # first, let's compute the CCFs between the 10s windows, for each pair
  tmp1 = wccfull[wccfull$triad==triads[i]&wccfull$typ=='LC'&wccfull$cond=='obs',]$wcc
  tmp2 = wccfull[wccfull$triad==triads[i]&wccfull$typ=='LR'&wccfull$cond=='obs',]$wcc
  x = ccf(tmp1,tmp2,lag.max=12,plot=F)$acf
  if (i==1) { # and we save the rows
    wccCcf = data.frame(triad=i,cond='obs',lag=-12:12,r=x)
  } else {
    wccCcf = rbind(wccCcf,data.frame(triad=i,cond='obs',lag=-12:12,r=x))
  }
  
  tmp1 = wccfull[wccfull$triad==triads[i]&wccfull$typ=='LC'&wccfull$cond=='obs',]$wcc
  tmp2 = wccfull[wccfull$triad==triads[i]&wccfull$typ=='CR'&wccfull$cond=='obs',]$wcc
  x = ccf(tmp1,tmp2,lag.max=12,plot=F)$acf
  wccCcf = rbind(wccCcf,data.frame(triad=i,cond='obs',lag=-12:12,r=x))
  
  tmp1 = wccfull[wccfull$triad==triads[i]&wccfull$typ=='LR'&wccfull$cond=='obs',]$wcc
  tmp2 = wccfull[wccfull$triad==triads[i]&wccfull$typ=='CR'&wccfull$cond=='obs',]$wcc
  x = ccf(tmp1,tmp2,lag.max=12,plot=F)$acf
  wccCcf = rbind(wccCcf,data.frame(triad=i,cond='obs',lag=-12:12,r=x))
  
}
par(mfrow=c(1,1))
# for plotting; assume low DF (# of triads, not pairs)
intactTriad = aggregate(r~lag,data=wccCcf,function(x) { c(m=mean(x),se=sd(x)/sqrt(35))})
wccCcfObs = wccCcf

# let's get the CCFs of the correlations: SURROGATE (vrt)
# precisely the same logic, but with the surrogates...
for (i in 1:35) {
  print(i)
  for (n in 1:34) { # we have n-1 surrogates for each triad
    tmp1 = wccfull[wccfull$triad==triads[i]&wccfull$typ=='LC'&wccfull$cond=='vrt'&wccfull$n==n,]$wcc
    tmp2 = wccfull[wccfull$triad==triads[i]&wccfull$typ=='LR'&wccfull$cond=='vrt'&wccfull$n==n,]$wcc
    x = ccf(tmp1,tmp2,lag.max=12,plot=F)$acf
    if (i==1) {
      wccCcf = data.frame(triad=i,cond='vrt',lag=-12:12,r=x)
    } else {
      wccCcf = rbind(wccCcf,data.frame(triad=i,cond='vrt',lag=-12:12,r=x))
    }
    
    tmp1 = wccfull[wccfull$triad==triads[i]&wccfull$typ=='LC'&wccfull$cond=='vrt'&wccfull$n==n,]$wcc
    tmp2 = wccfull[wccfull$triad==triads[i]&wccfull$typ=='CR'&wccfull$cond=='vrt'&wccfull$n==n,]$wcc
    x = ccf(tmp1,tmp2,lag.max=12,plot=F)$acf
    wccCcf = rbind(wccCcf,data.frame(triad=i,cond='vrt',lag=-12:12,r=x))
    
    tmp1 = wccfull[wccfull$triad==triads[i]&wccfull$typ=='LR'&wccfull$cond=='vrt'&wccfull$n==n,]$wcc
    tmp2 = wccfull[wccfull$triad==triads[i]&wccfull$typ=='CR'&wccfull$cond=='vrt'&wccfull$n==n,]$wcc
    x = ccf(tmp1,tmp2,lag.max=12,plot=F)$acf
    wccCcf = rbind(wccCcf,data.frame(triad=i,cond='vrt',lag=-12:12,r=x)) 
  }
}

par(mfrow=c(1,1))

# for plotting; assume low DF (# of triads, not surrogates)
surrogateTriad = aggregate(r~lag,data=wccCcf,function(x) { c(m=mean(x),se=sd(x)/sqrt(35))})
# put them together (observed / surrogate)
wccCcf = rbind(wccCcfObs,wccCcf)

# statistical test at lag 0
lmo = lmer(r~cond+(1+cond|triad),data=wccCcf[wccCcf$lag==0,])
pander(print_stats(lmo))

# let's plot 'em now
source('plotCCF.R')

#
# exploratory regression analysis
# 
source('combineMansonData.R') # integrate with manson task/individual diff'ces data
# note it stores only one of LC/CL, etc., averaging across positions
# NB: makes mansonData

# again, treat as large repeated measures with a single intercept
lmo = lmer(wcc~as.factor(sex)+p1.income.zip+
             p1.primary.psycho+p2.facial+
             cultural.style+
             language.style.match+
             common.ground+
             p2.interrupts.p1+
             p1.pd.toward.p2+
             p1.rates.p2.warmth+
             p1.rates.p2.competence+
             laughDat$Total.laughs+
             laughDat$X..SHARED+(1|triad),data=wccMansonData)
print_stats(lmo)

# let's make sure these square away even with a full maximal model
lmo = lmer(wcc~cultural.style+language.style.match+laughDat$X..SHARED+(1+cultural.style+language.style.match+laughDat$X..SHARED|triad),data=wccMansonData)
print_stats(lmo)

# if there's compensation, then there should be significant prediction by interaction
ccfmaxC = scale(mansonData$ccfmax) # scale all variables to center, and normalize
lingstyleC = scale(mansonData$language.style.match)
lmo = glmer(p1.pd.toward.p2~ccfmaxC*lingstyleC+(1|triad),data=mansonData,family=binomial)
summary(lmo)

# what's the interaction; let's plot the apparently strongest interaction
pdf('figures/figure-5.pdf',height=5,width=9)
par(mfrow=c(1,2))

ixes = which(mansonData$p1.pd.toward.p2==1)
plot((lingstyleC[ixes]),
     (ccfmaxC[ixes]),xlab='z LSM',ylab='z max[CCF]',
     type='p',col=rgb(0,0,0,.5),pch=1,cex=1,main='Cooperate',ylim=c(-2,4))
abline(lm(ccfmaxC[ixes]~lingstyleC[ixes]),lwd=3,col='black',lty='dashed')

ixes = which(mansonData$p1.pd.toward.p2==0)
plot((lingstyleC[ixes]),
     (ccfmaxC[ixes]),xlab='z LSM',ylab='z max[CCF]',
     type='p',col=rgb(0,0,0,.5),pch=1,cex=1,main='Defect',ylim=c(-2,4))
abline(lm(ccfmaxC[ixes]~lingstyleC[ixes]),lwd=3,col='black',lty='dashed')

dev.off()

