#
# this segment of code loops through all motion signals from the dats
# entries and computes a series of DVs on the motion, including
# windowed cross correlation (saved in ccfres)
# note: flst is defined in bodyExtract.R, and dats is from the data file that script saves
#

wsz = 300
wshft = 300
ccfres = c(); wccres = c(); wccfull = c()
for (fl in flst) {
  print(paste('Processing triad in:',fl))
  otherfs = sample(flst[flst!=fl],length(flst)-1) # get list of others, for baseline
  tgt = dats[dats$triad==fl,] # get data from this triad
  tgt = tgt[200:dim(tgt)[1],] 
  # the low-pass filter always induces a rise in amplitude in the first part of the time series; this could produce spurious correlation, so we remove
  
  #
  # now get the windowed cross correlation and overall cross correlation function
  # between each pair of individuals in the triad
  #
  ccfObsLC = as.numeric(ccf(tgt$filterLeft,tgt$filterMid,plot=F,lag.max=100)$acf) 
  ccWinLC = windowedcc(tgt$filterLeft,tgt$filterMid,wsz,wshft) 
  
  ccfObsCR = as.numeric(ccf(tgt$filterMid,tgt$filterRight,plot=F,lag.max=100)$acf) 
  ccWinCR = windowedcc(tgt$filterMid,tgt$filterRight,wsz,wshft)
  
  ccfObsLR = as.numeric(ccf(tgt$filterLeft,tgt$filterRight,plot=F,lag.max=100)$acf) 
  ccWinLR = windowedcc(tgt$filterLeft,tgt$filterRight,wsz,wshft) 
  
  # let's aggregate results
  ccfres = rbind(ccfres,data.frame(triad=fl,ccf=ccfObsLC,typ='LC',cond='obs',lag=-100:100))
  ccfres = rbind(ccfres,data.frame(triad=fl,ccf=ccfObsCR,typ='CR',cond='obs',lag=-100:100))
  ccfres = rbind(ccfres,data.frame(triad=fl,ccf=ccfObsLR,typ='LR',cond='obs',lag=-100:100))
  
  # contains max corr, max lag location, etc.
  wccres = rbind(wccres,data.frame(triad=fl,cond='obs',typ='LC',
                             winmin=min(ccWinLC),winmax=max(ccWinLC),
                             ccfmax=max(ccfObsLC),winsd=sd(ccWinLC),
                             max.loc=which.max(ccfObsLC)))
  wccres = rbind(wccres,data.frame(triad=fl,cond='obs',typ='CR',
                                   winmin=min(ccWinCR),winmax=max(ccWinCR),
                                   ccfmax=max(ccfObsCR),winsd=sd(ccWinCR),
                                   max.loc=which.max(ccfObsCR)))
  wccres = rbind(wccres,data.frame(triad=fl,cond='obs',
                                   typ='LR',winmin=min(ccWinLR),
                                   winmax=max(ccWinLR),ccfmax=max(ccfObsLR),
                                   winsd=sd(ccWinLR),
                                   max.loc=which.max(ccfObsLR)))
  
  # contains raw windowed correlation, row wise, across all
  # the n variable is used to track n surrogates, below; here, n is just 1 (1 observed)
  wccfull = rbind(wccfull,data.frame(triad=fl,cond='obs',typ='LC',wcc=ccWinLC,t=1:length(ccWinLC),n=1))
  wccfull = rbind(wccfull,data.frame(triad=fl,cond='obs',typ='CR',wcc=ccWinCR,t=1:length(ccWinCR),n=1))
  wccfull = rbind(wccfull,data.frame(triad=fl,cond='obs',typ='LR',wcc=ccWinLR,t=1:length(ccWinLR),n=1))
  
  #
  # now time to build a surrogate baseline by building cc and wcc from 
  # pseudo pairs across the triads
  # the process is precisely the same, but we use comp / tgt instead
  # for completeness, we obtain all such surrogates 
  #  
  for (otherfl in otherfs) {    
    n=which(otherfs==otherfl)
    comp = dats[dats$triad==otherfl,]
    comp = comp[200:dim(comp)[1],]
    
    ccfObsLC = as.numeric(ccf(comp$filterLeft,tgt$filterMid,plot=F,lag.max=100)$acf) ####### ts types
    ccWinLC = windowedcc(comp$filterLeft,tgt$filterMid,wsz,wshft) ####### ts types
    ccfObsCR = as.numeric(ccf(comp$filterMid,tgt$filterRight,plot=F,lag.max=100)$acf) ####### ts types
    ccWinCR = windowedcc(comp$filterMid,tgt$filterRight,wsz,wshft) ####### ts types
    ccfObsLR = as.numeric(ccf(comp$filterLeft,tgt$filterRight,plot=F,lag.max=100)$acf) ####### ts types
    ccWinLR = windowedcc(comp$filterLeft,tgt$filterRight,wsz,wshft) ####### ts types
    
    wccres = rbind(wccres,data.frame(triad=fl,cond='vrt',typ='LC',
                               winmin=min(ccWinLC),winmax=max(ccWinLC),
                               ccfmax=max(ccfObsLC),winsd=sd(ccWinLC),max.loc=which.max(ccfObsLC)))
    wccres = rbind(wccres,data.frame(triad=fl,cond='vrt',typ='CR',
                               winmin=min(ccWinCR),winmax=max(ccWinCR),
                               ccfmax=max(ccfObsCR),winsd=sd(ccWinCR),max.loc=which.max(ccfObsCR)))
    wccres = rbind(wccres,data.frame(triad=fl,cond='vrt',typ='LR',
                               winmin=min(ccWinLR),winmax=max(ccWinLR),
                               ccfmax=max(ccfObsLR),winsd=sd(ccWinLR),max.loc=which.max(ccfObsLR)))
    
    ccfres = rbind(ccfres,data.frame(triad=fl,ccf=ccfObsLC,typ='LC',cond='vrt',lag=-100:100))    
    ccfres = rbind(ccfres,data.frame(triad=fl,ccf=ccfObsCR,typ='CR',cond='vrt',lag=-100:100))
    ccfres = rbind(ccfres,data.frame(triad=fl,ccf=ccfObsLR,typ='LR',cond='vrt',lag=-100:100))
    
    wccfull = rbind(wccfull,data.frame(triad=fl,cond='vrt',typ='LC',wcc=ccWinLC,t=1:length(ccWinLC),n=n))
    wccfull = rbind(wccfull,data.frame(triad=fl,cond='vrt',typ='CR',wcc=ccWinCR,t=1:length(ccWinCR),n=n))
    wccfull = rbind(wccfull,data.frame(triad=fl,cond='vrt',typ='LR',wcc=ccWinLR,t=1:length(ccWinLR),n=n))
    
  }
  
  plotCcfs(fl,wccres,ccfres)
  
}

save(ccfres,file='ccfresData.RData')
save(wccres,file='wccData.RData')
save(wccfull,file='wccfullData.RData')




