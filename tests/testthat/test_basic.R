library(Hmisc)

set.seed(2)
data<-read.csv("~/Dropbox/Berkeley_Projects/School Projects/COMPSCI294_F2017_Project/data/data.csv", row.names = 1)

fit<-initEst(data, freqW = 2,freqA = 2,freqY = 2)

#Intervention at time point 12
mcEst_int <- mcEst(fit, start=1, node="W", t=5, Anode=3, intervention=1, lag=0, MC=1000)
mcEst_int$estimate

#No intervention
mcEst_noint <- mcEst(fit, start=1, node="W", t=5, Anode=3, intervention=NULL, lag=0, MC=1000)
mcEst_noint$estimate

#Calculate the clever covariate
cleverCovariate<-cleverCov(fit,t=5,Anode=3,MC=10)




