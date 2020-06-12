
##
# function extracts the time series from the optic flow analyzer, filters them, plots some
# first: set working directory to this file's location
##

source('functions.R')
library(gdata)
library(signal)

flst = read.table('raw_flow_data/videoList.txt',stringsAsFactors=F)$V1

for (fl in flst) { # let's load up all the data into a data frame
  print(fl)
  if (fl == flst[1]) {
    dats = getData(paste('raw_flow_data/',fl,sep=''),fl,plotit=T)
    dats$triad = fl
  } else {
    dats_temp = getData(paste('raw_flow_data/',fl,sep=''),fl,plotit=T)  
    dats_temp$triad = fl
    dats = rbind(dats,dats_temp)
  }
}
colnames(dats) = list('filterLeft','filterMid','filterRight','triad')

save(file='loadedSignals.RData',dats)