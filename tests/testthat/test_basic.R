library(Hmisc)

set.seed(2)
data<-read.csv("~/Dropbox/Berkeley_Projects/School Projects/COMPSCI294_F2017_Project/data/data.csv", row.names = 1)
names<-row.names(data)

#Make it a bit easier:
data<-data.frame(data=data[1:303,])
row.names(data)<-names[1:303]

#Get initial estimates of the process:
fit<-initEst(data, freqW = 2,freqA = 2,freqY = 2)

#Calculate the clever covariate with default settings (B=N=MC=100)
#cleverCovariate<-cleverCov(fit,t=5,Anode=3, intervention=1)

#Get estimates:
est<-mainTMLE(fit, t=5, Anode=5)





