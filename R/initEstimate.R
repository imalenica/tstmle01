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
#' }         
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
    fitW<-glm(formula = data ~ ., family = quasibinomial(link = "logit"), data = W)
    fitA<-glm(formula = data ~ ., family = quasibinomial(link = "logit"), data = A)
    fitY<-glm(formula = data ~ ., family = quasibinomial(link = "logit"), data = Y)
    
    fit<-list(W=fitW,A=fitA,Y=fitY,freqW=freqW,freqA=freqA,freqY=freqY,data=data,lag_data=data_est)
    
  }
  
  return(fit)
}
