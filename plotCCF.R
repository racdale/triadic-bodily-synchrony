# let's plot the ccf function for intact / surrogate
plotfl = paste('figures/figure-4.pdf',sep='')
pdf(file=plotfl,height=3.5,width=10)
par(mfrow=c(1,3),mar=c(4,4,2,2))

intact = aggregate(ccf~lag,
                   data=ccfres[ccfres$cond=='obs',],
                   function(x) { c(m=mean(x),se=sd(x)/sqrt(35))})

plot(intact$lag*.033,intact$ccf[,1],type='l', # 33ms sample rate
     xlab='Relative lag (s)',main='Dyad cross-correlation function',
     ylab='Correlation coefficient (r)',
     ylim=c(-.03,.08),lwd=3,col='green')
points(intact$lag*.033,intact$ccf[,1]+intact$ccf[,2],type='l',col='green')
points(intact$lag*.033,intact$ccf[,1]-intact$ccf[,2],type='l',col='green')

surrogate = aggregate(ccf~lag, 
                      data=ccfres[ccfres$cond=='vrt',],
                      function(x) { c(m=mean(x),se=sd(x)/sqrt(35))})
points(surrogate$lag*.033,surrogate$ccf[,1],type='l',col='red',lwd=3)
points(surrogate$lag*.033,surrogate$ccf[,1]+surrogate$ccf[,2],type='l',col='red')
points(surrogate$lag*.033,surrogate$ccf[,1]-surrogate$ccf[,2],type='l',col='red')

lagLocs = wccres[wccres$cond=='obs',]
hist(lagLocs$max.loc*.033,5,main='Maximum lag location distribution',
     xlab='Maximum lag location (s)',xlim=c(-6,6),
     ylab='Number of dyads')

plot(intactTriad$lag,intactTriad$r[,1],type='l',
     xlab='Lag (10s window)',main='Triad cross-correlation function',
     ylab='Correlation coefficient (r)',
     ylim=c(-.04,.2),lwd=3,col='green')
points(intactTriad$lag,intactTriad$r[,1]+intactTriad$r[,2],type='l',col='green')
points(intactTriad$lag,intactTriad$r[,1]-intactTriad$r[,2],type='l',col='green')

points(surrogateTriad$lag,surrogateTriad$r[,1],type='l',col='red',lwd=3)
points(surrogateTriad$lag,surrogateTriad$r[,1]+surrogateTriad$r[,2],type='l',col='red')
points(surrogateTriad$lag,surrogateTriad$r[,1]-surrogateTriad$r[,2],type='l',col='red')

dev.off()

