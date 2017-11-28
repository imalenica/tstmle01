###############################
#Simulate binary time series.
###############################

library(simcausal)
options(simcausal.verbose=FALSE)

t.end <- 10

#No missing values, keep it very simple
D <- DAG.empty()
D <- D +
  
  #Initialize first points with probability 0.5?
  node("W", t = 0, distr = "rbern",prob = 0.5) +
  node("A", t = 0, distr = "rbern",prob = 0.5) +
  node("Y", t = 0, distr = "rbern",prob = 0.5)
  
# Define very simple dependence- this is technically not even order(1)?
# W(t) should depend on Y(t-1) and A(t-1)
# A(t) should depend on W(t) and Y(t-1)
# Y(t) should depend on A(t) and W(t)

D <- D +
  node("W", t = 1:t.end, distr = "rbern",
       prob = plogis(-5 + 4*A[t-1] + 3*Y[t-1])) +
  node("A", t = 1:t.end, distr = "rbern",
       prob = plogis(0.2 + 0.7*W[t] - 0.2*Y[t-1])) +
  node("Y", t = 1:t.end, distr = "rbern",
       prob = plogis(-1.2 + 0.7*W[t]+ 1/3*A[t]))

lDAG <- set.DAG(D)

# plot the DAG for different tmax
plotDAG(lDAG, tmax = 3, xjitter = 0.3, yjitter = 0.03,
        edge_attrs = list(width = 0.5, arrow.width = 0.4, arrow.size = 0.8),
        vertex_attrs = list(size = 12, label.cex = 0.8))

plotDAG(lDAG, tmax = 10, xjitter = 0.3, yjitter = 0.03,
        edge_attrs = list(width = 0.5, arrow.width = 0.4, arrow.size = 0.8),
        vertex_attrs = list(size = 12, label.cex = 0.8))

# take a draw of 1 from this DGP
data1 <- sim(DAG = lDAG, n = 1, rndseed = 123)

###############################################################
#Include dependence on past self to make it actually order(1)
###############################################################

D <- DAG.empty()
D <- D +
  
  #Initialize first points with probability 0.5?
  node("W", t = 0, distr = "rbern",prob = 0.5) +
  node("A", t = 0, distr = "rbern",prob = 0.5) +
  node("Y", t = 0, distr = "rbern",prob = 0.5)

# W(t) should depend on Y(t-1), A(t-1), W(t-1)
# A(t) should depend on W(t), Y(t-1), A(t-1)
# Y(t) should depend on A(t), W(t), Y(t-1)

D <- D +
  node("W", t = 1:t.end, distr = "rbern",
       prob = plogis(-5 + 4*A[t-1] + 3*Y[t-1] + 0.2*W[t-1])) +
  node("A", t = 1:t.end, distr = "rbern",
       prob = plogis(0.2 + 0.7*W[t] - 0.2*Y[t-1] + 0.3*A[t-1])) +
  node("Y", t = 1:t.end, distr = "rbern",
       prob = plogis(-1.2 + 0.7*W[t]+ 1/3*A[t] - 0.5*Y[t-1]))

lDAG <- set.DAG(D)

# plot the DAG for different tmax
plotDAG(lDAG, tmax = 3, xjitter = 0.3, yjitter = 0.03,
        edge_attrs = list(width = 0.5, arrow.width = 0.4, arrow.size = 0.8),
        vertex_attrs = list(size = 12, label.cex = 0.8))

plotDAG(lDAG, tmax = 10, xjitter = 0.3, yjitter = 0.03,
        edge_attrs = list(width = 0.5, arrow.width = 0.4, arrow.size = 0.8),
        vertex_attrs = list(size = 12, label.cex = 0.8))

# take a draw of 1 from this DGP
data2<- sim(DAG = lDAG, n = 1, rndseed = 123)
