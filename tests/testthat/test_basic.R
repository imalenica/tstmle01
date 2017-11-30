
fit<-initEst(data,freqW = 2,freqA = 2,freqY = 2)
mcEst <- mcEst(fit, t=20, Anode=12, intervention=1, MC=100)
mcEst <- mcEst(fit, t=20, Anode=12, intervention=NULL, MC=100)
