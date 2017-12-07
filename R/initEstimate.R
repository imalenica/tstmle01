#' Initial Estimate
#'
#' This function generates coefficients of previous variables based on a single time series observation.
#' Later on, this should include more complex relationships instead of just past variables and 
#' better estimation methods.
#'
#' @param data data.frame object containg the time series with relevant time ordering. 
#' @param freqW A numeric specifying the Markov order for W nodes. 
#' @param freqA A numeric specifying the Markov order for A nodes. 
#' @param freqY A numeric specifying the Markov order for Y nodes. 
#'
#' @return An object of class \code{tstmle}.
#' \describe{
#' \item{fitW}{Fit object for W part of the likelihood.}
#' \item{fitA}{Fit object for A part of the likelihood.}
#' \item{fitY}{Fit object for Y part of the likelihood.}  
#' \item{freqW}{Input Markov order for W.}
#' \item{freqA}{Input Markov order for A.}
#' \item{freqY}{Input Markov order for Y.}
#' \item{data}{Input data.}
#' \item{lag_data}{Data with lagged components.}
#'         
#'                         
#' @importFrom Hmisc Lag
#' 
#'
#' @export
#'

initEst <- function(data, freqW=NULL, freqA=NULL, freqY=NULL) {
  
  if(all.equal(freqA,freqW,freqY)){
    
    #Lag past 
    res<-lapply(1:freqW, function(x) {Lag(data[,1], x)})
    res <- data.frame(matrix(unlist(res), nrow=length(res[[1]])), stringsAsFactors = FALSE)
    
    data_est_full<-cbind.data.frame(data=data,res)

    #Drop time 0 for estimation
    cc<-complete.cases(data_est_full)
    data_est<-data_est_full[cc,]
    data_est<-data.frame(data_est[-1,])
    
    #Separate by process:
    W<-data_est[grep('W', row.names(data_est), value=TRUE),]
    A<-data_est[grep('A', row.names(data_est), value=TRUE),]
    Y<-data_est[grep('Y', row.names(data_est), value=TRUE),]
    
    #Estimate the process coefficients. This should allow better strategies later on
    fitW<-glm(formula = data ~ ., family = binomial(link = "logit"), data = W)
    fitA<-glm(formula = data ~ ., family = binomial(link = "logit"), data = A)
    fitY<-glm(formula = data ~ ., family = binomial(link = "logit"), data = Y)
    
    fit<-list(W=fitW,A=fitA,Y=fitY,freqW=freqW,freqA=freqA,freqY=freqY,data=data,lag_data=data_est)
    
  }else if(!all.equal(freqA,freqW,freqY)){
    
    resW <- lapply(1:freqW, function(x) {Lag(data[,1], x)})
    resW <- data.frame(matrix(unlist(resW), nrow=length(resW[[1]])), stringsAsFactors = FALSE)
    data_estW<-cbind.data.frame(data=data,resW)
    
    #Drop time 0 for estimation
    cc<-complete.cases(data_estW)
    data_estW<-data_estW[cc,]
    data_estW<-data.frame(data_estW[-1,])
    
    resA <- lapply(1:freqA, function(x) {Lag(data[,1], x)})
    resA <- data.frame(matrix(unlist(resA), nrow=length(resA[[1]])), stringsAsFactors = FALSE)
    data_estA<-cbind.data.frame(data=data,resA)
    
    #Drop time 0 for estimation
    cc<-complete.cases(data_estA)
    data_estA<-data_estA[cc,]
    data_estA<-data.frame(data_estA[-1,])
    
    resY <- lapply(1:freqY, function(x) {Lag(data[,1], x)})
    resY <- data.frame(matrix(unlist(resY), nrow=length(resY[[1]])), stringsAsFactors = FALSE)
    data_estY<-cbind.data.frame(data=data,resY)
    
    #Drop time 0 for estimation
    cc<-complete.cases(data_estY)
    data_estY<-data_estY[cc,]
    data_estY<-data.frame(data_estY[-1,])
    
    W<-data_est[grep('W', row.names(data_estW), value=TRUE),]
    A<-data_est[grep('A', row.names(data_estA), value=TRUE),]
    Y<-data_est[grep('Y', row.names(data_estY), value=TRUE),]
    
    #Estimate the process coefficients. This should allow better strategies later on
    fitW<-glm(formula = data ~ ., family = binomial(link = "logit"), data = W)
    fitA<-glm(formula = data ~ ., family = binomial(link = "logit"), data = A)
    fitY<-glm(formula = data ~ ., family = binomial(link = "logit"), data = Y)
    
    fit<-list(W=fitW,A=fitA,Y=fitY,freqW=freqW,freqA=freqA,freqY=freqY,data=data)
    
  }else{
    
    fit<-forecast::auto.arima(data)
    
  }
  
  return(fit)
}