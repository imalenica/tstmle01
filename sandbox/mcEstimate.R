#' Compute the Monte Carlo estimate
#'
#' This function computes the Monte Carlo estimate of the expected value of the parameter of interest
#' under specified intervention and under the estimated mechanism. 
#' 
#' @param fit \code{fit} object obtained by \code{initEst}. 
#' @param t Outcome time point of interest 
#' @param Anode Intervention node.
#' @param intervention Specify g^*, of P(A|past).   
#'
#' @return An object of class \code{tstmle}.
#' \describe{
#' \item{fit}{Fit objects for each part of the likelihood.}
#'            
#'
#'
#' @export
#'

mcEst <- function(fit, t, Anode, intervention, MC) {

  data<-fit$data
  
  #How many in a batch:
  step<-length(grep('_0', row.names(data), value=TRUE))
  
  data_lag<-fit$lag_data
  data_lag<-data.frame(data_lag[,-1])
  data_lag<-test[1:(t*step),]
  estNames<-row.names(data_lag)
  
  #get actual index of Anode (based on data_lag)
  Anode<-Anode*step-1
  
  #Later on will need to include more Ws
  #This should be parallelized 
  outcome<-matrix(nrow = MC, ncol = 1)
  data_lag_origin<-data_lag
  
  for(B in 1:MC){
    
    data_lag<-data_lag_origin
    
    for(i in seq(1,nrow(data_lag),step)){
      
      newW<-rbinom(1,1,plogis(predict(fit$W, data_lag[i,], type="response")))
      data_lag[(i+1),1]<-newW
      if(i != 1){
        data_lag[(i+1),2]<-newY
      }
      
      if(is.null(intervention)){
        newA<-rbinom(1,1,plogis(predict(fit$A, data_lag[i+1,], type="response")))
        data_lag[(i+2),2]<-newW
        data_lag[(i+2),1]<-newA
      }else if((i+1) == Anode){
        newA<-rbinom(1,1,intervention)
        data_lag[(i+2),2]<-newW
        data_lag[(i+2),1]<-newA
      }
      
      newY<-rbinom(1,1,plogis(predict(fit$Y, data_lag[i+2,], type="response")))
      data_lag[(i+3),2]<-newA
      data_lag[(i+3),1]<-newY
      
    }
    
    outcome[B,]<-newY
    
  }
  
  return(list(estimate=round(mean(outcome)), outcome=outcome, intervention=intervention,MC=MC, t=t, Anode=Anode))

 
}

