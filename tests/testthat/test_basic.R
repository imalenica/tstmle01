library(Hmisc)
library(doParallel)
library(parallel)
library(foreach)

set.seed(2)
data("ts_samp_data")
data<-ts_samp_data
names<-row.names(data)

# Make it very easier:
data <- data.frame(data = data[1:303, ])
row.names(data) <- names[1:303]

# Get initial estimates of the process:
fit <- initEst(data, t=5, freqW = 2, freqA = 2, freqY = 2)

# Get estimates:
est1 <- mainTMLE(fit, t = 5, Anode = 3, intervention = 1, B = 100, N = 100, MC = 100, maxIter = 20)
est2 <- mainTMLE(fit, t = 5, Anode = 3, B = 100, N = 100, MC = 50, maxIter = 20)
