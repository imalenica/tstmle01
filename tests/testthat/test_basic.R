library(Hmisc)

set.seed(2)
data <- read.csv("~/Dropbox/Berkeley_Projects/Software/tstmle01/data/data.csv", row.names = 1)
names <- row.names(data)

# Make it very easy:
data <- data.frame(data = data[1:33, ])
row.names(data) <- names[1:33]

# Get initial estimates of the process:
fit <- initEst(data, freqW = 2, freqA = 2, freqY = 2)

# Get estimates:
est1 <- mainTMLE(fit, t = 5, Anode = 3, intervention = 1, B = 50, N = 50, MC = 50, maxIter = 20)
est2 <- mainTMLE(fit, t = 5, Anode = 3, B = 50, N = 50, MC = 50, maxIter = 20)
