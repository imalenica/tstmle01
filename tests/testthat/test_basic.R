library(Hmisc)

set.seed(2)
data<-read.csv("~/Dropbox/Berkeley_Projects/Software/tstmle/data/data.csv", row.names = 1)

fit<-initEst(data, freqW = 2,freqA = 2,freqY = 2)

#Intervention at time point 12
mcEst_int <- mcEst(fit, start=1, node="W", t=20, Anode=12, intervention=1, MC=1000)
mcEst_int$estimate

#No intervention
mcEst_noint <- mcEst(fit, start=1, node="W", t=20, Anode=12, intervention=NULL, MC=1000)
mcEst_noint$estimate






