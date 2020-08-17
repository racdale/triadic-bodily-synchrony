# first: set working directory to here

mansonDataAll = read.csv('Manson_Convo_PD_data.csv',header=T,skip=1)
colnames(mansonDataAll) = list('order','triad','code','chairs','p1.id','sex','p1.income.zip','p1.primary.psycho','p1.secondary.psycho','p2.facial','cultural.style','language.style.match','common.ground','p2.interrupts.p1','p1.pd.toward.p2','p2.pd.toward.p1','p1.rates.p2.warmth','p1.rates.p2.competence','p2.rates.p1.warmth','p2.rates.p1.competence')

mansonData = mansonDataAll 
# in case we want to focus analysis on one chair ordering
# mansonData = mansonDataAll[mansonDataAll$chairs=='LC'|mansonDataAll$chairs=='CR'|mansonDataAll$chairs=='LR',] 

laughs = read.table('laughter_data.txt',header=T,sep='\t') # laughter data from bryant, several available columns, pre-selected 2 obvious ones

mansonData$winmin = -1 # initialize the columns we will populate with the body DVs
mansonData$winmax = -1
mansonData$ccfmax = -1
mansonData$winsd = -1
mansonData$max.loc = -1

mansonData$chairs = as.character(mansonData$chairs) # chair ordering in manson
wccres$typ = as.character(wccres$typ) # typ = chair ordering in body data

wccresNew = wccres[wccres$cond=='obs',] # store observed data (without surrogate)

reverseString = function(x) { # for checking measures from reverse-chair rows in manson data
  return(paste(substr(x,2,2),substr(x,1,1),sep='')) 
} # e.g. subjects specked as LR in all data; however we want RL in manson, since we want to know what subject R thinks about warmth, PD, etc.

mansonLaughs = c() # build this separately for inclusion later
for (i in 1:dim(wccresNew)[1]) { # integrating Manson data with body correlations (windowed correlation = wcc)
  ix=regexpr('_',wccresNew[i,]$triad)[1]-1
  triad = as.numeric(substr(wccresNew[i,]$triad,2,ix)) # get triad # from body movement file, e.g., T8_...
  mansonLaughs = rbind(mansonLaughs,laughs[laughs$CONV==triad,2:3],laughs[laughs$CONV==triad,2:3]) # get the laughs
  # get measures from both orderings in manson data (LR / RL)
  mansonData[mansonData$chairs==wccresNew[i,]$typ & mansonData$triad==triad,21:25] = wccresNew[i,4:8]
  mansonData[mansonData$chairs==reverseString(wccresNew[i,]$typ) & mansonData$triad==triad,21:25] = wccresNew[i,4:8]
}

#
# let's retrieve the windowed correlation scores and treat this as a repeated
# measures setup... multiple observations in 10s segments
# increases power for the exploratory analysis as described in the main paper
#
wccMansonData = wccfull[wccfull$cond=='obs',]
desiredCols = c('p1.id','sex','p1.income.zip','p1.primary.psycho','p2.facial',
                'cultural.style','language.style.match','common.ground',
                'p2.interrupts.p1','p1.pd.toward.p2','p1.rates.p2.warmth',
                'p1.rates.p2.competence')
lapply(desiredCols,function(x) {
  thisExpr = paste('wccMansonData$',unlist(x),'<<- -99',sep='')
  eval(parse(text=thisExpr))
})

l = nrow(wccMansonData)
wccMansonDataRev = wccMansonData # so we can get the reverse-chair scores
laughDat = c()

#
# a slow and lame loop... but we just do it once and we're done
# *NB: wccMansonData = windowed correlations combined with Manson et al. covariates
#
for (i in 1:l) {
  print(i/l)
  ix=regexpr('_',wccMansonData[i,]$triad)[1]-1
  triad = as.numeric(substr(wccMansonData[i,]$triad,2,ix))
  # we do laughter as a separate vector since it came from another analysis of the videos
  # by aligning the vector, we can just slip it into the covariate list in the model
  laughDat = rbind(laughDat,laughs[laughs$CONV==triad,2:3]) 
  typ = wccMansonData[i,]$typ
  revTyp = reverseString(wccMansonData[i,]$typ)
  mansonVect = mansonData[mansonData$chairs==typ&mansonData$triad==triad,]
  wccMansonData[i,7:18] = subset(mansonVect,select=desiredCols)
  mansonVect = mansonData[mansonData$chairs==revTyp&mansonData$triad==triad,]
  wccMansonData[i,7:18] = (wccMansonData[i,7:18] + subset(mansonVect,select=desiredCols))/2
}

#
# let's check for Table 1 in paper... make sure we get the right N for SD:
#
tableData = mansonData[mansonData$chairs%in%c('LC','LR','CR'),]
summary(tableData)
sd(laughs$X..SHARED*100) # checking to confirm laughter rows in Table 1


