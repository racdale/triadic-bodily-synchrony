
#
#  print out statistics from lmer
#
print_stats = function(lmo) {
  coefs = data.frame(summary(lmo)$coefficient)
  coefs$p = 2*(1-pnorm(abs(coefs$t.value)))
  return(coefs)
}

#
# calculate windowed cross correlation (a-la Boker et al.)
# input x time series, y time series, window size, and shift size
# outputs sequence of lag-0 correlations at that window size
#
windowedcc = function(x,y,wsz,shft) {
  x = x[1:(min(length(x),length(y)))]
  y = y[1:(min(length(x),length(y)))]
  cors = c()
  ixes = seq(from=1,to=(length(x)-wsz+1),by=shft)
  for (i in ixes) {
    xt = x[i:(i+wsz-1)]
    yt = y[i:(i+wsz-1)]
    cors = c(cors,cor(xt,yt))
  }
  return(cors)
}

#
# get the optic flow from the vatikiotis-bateson output
# generates low-pass filtered version of time series located
# in the xls file (x, y, absolute pixel change)
#
getData = function(flpath,fl,plotit=F) {
  a = read.xls(flpath)
  backup_a = a
  header = a[1,]
  a = a[2:dim(a)[1],]
  newColNames = c()
  for (i in 1:length(header)) { # clean column names
    newColNames = c(newColNames,paste(colnames(a)[i],as.character(header[,i]),sep=''))
  }
  colnames(a) = newColNames
  for (i in 1:dim(a)[2]) { 
    a[,i] = as.numeric(as.character(a[,i]))
  }
  bf = butter(8,.05,'low')
  #fz=freqz(bf,Fs=30)
  filterMid = filter(bf,a$middlemag)
  
  bf = butter(8,.05,'low')
  filterLeft = filter(bf,a$leftmag)
  
  bf = butter(8,.05,'low')
  filterRight = filter(bf,a$rightmag)

  if (plotit) {
    plotfl = paste('figures/sample_time_series/',fl,'.pdf',sep='')
    pdf(file=plotfl,height=5,width=5)
    par(mfrow=c(3,1),mar=c(4,4,2,2))
    plot(a$Time..s.NA[1:2000],a$leftmag[1:2000],
         col='grey',type='p',xlab='Time (s)',ylab='Magnitude change',main='Left member')
    # alignment requires a 1s shift on filter...
    points(a$Time..s.NA[1:2000],filterLeft[31:2030],col='black',type='l')
    plot(a$Time..s.NA[1:2000],a$middlemag[1:2000],
         col='grey',type='p',xlab='Time (s)',ylab='Magnitude change',main='Center member')
    points(a$Time..s.NA[1:2000],filterMid[31:2030],col='black',type='l')
    plot(a$Time..s.NA[1:2000],a$rightmag[1:2000],
         col='grey',type='p',xlab='Time (s)',ylab='Magnitude change',main='Right member')
    points(a$Time..s.NA[1:2000],filterRight[31:2030],col='black',type='l')        
    dev.off()
  }
  
  return(data.frame(as.numeric(filterLeft),as.numeric(filterMid),as.numeric(filterRight)))
}


plotCcfs = function(fl,wccres,ccfres) {
  subwcc = wccfull[wccfull$triad==fl&wccfull$cond=='obs',] # windowed correlation, over 300 samples
  subccf = ccfres[ccfres$triad==fl&ccfres$cond=='obs',] # cross correlation function across whole time series
  plotfl = paste('figures/dvs/',fl,'.pdf',sep='')
  pdf(file=plotfl,height=5,width=5)
  
  par(mfrow=c(2,2),mar=c(4,4,2,2))
  dt = subwcc[subwcc$typ=='LC',]
  plot(dt$t/60,dt$wcc,
       ylab='Correlation (r)',main='Left / Center Members',
       type='l',xlab='Time window (min)')
  points(which.max(dt$wcc)/60,max(dt$wcc),pch=15,cex=2)
  points(which.min(dt$wcc)/60,min(dt$wcc),pch=15,cex=2)
  dt = subccf[subccf$typ=='LC'&subccf$cond=='obs',]
  plot(dt$lag*.033,dt$ccf,
       ylab='Correlation (r)',main='Left / Center Members',
       type='l',xlab='Relative lag (s)')
  points((which.max(dt$ccf)-100)*.033,max(dt$ccf),pch=15,cex=2)
  dt = subwcc[subwcc$typ=='CR',]
  plot(dt$t/60,dt$wcc,
       ylab='Correlation (r)',main='Center / Right Members',
       type='l',xlab='Time window (10 s)')
  points(which.max(dt$wcc)/60,max(dt$wcc),pch=15,cex=2)
  points(which.min(dt$wcc)/60,min(dt$wcc),pch=15,cex=2)
  dt = subccf[subccf$typ=='CR'&subccf$cond=='obs',]
  plot(dt$lag*.033,dt$ccf,
       ylab='Correlation (r)',main='Center / Right Members',
       xlab='Relative lag (s)',type='l')
  points((which.max(dt$ccf)-100)*.033,max(dt$ccf),pch=15,cex=2)
  dev.off()
    
}

