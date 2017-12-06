#' Compute the Monte Carlo estimate
#'
#' This function computes the Monte Carlo estimate of the expected value of the parameter of interest
#' under specified intervention and under the estimated mechanism. 
#' 
#' @param fit \code{fit} object obtained by \code{initEst}. 
#' @param start Start generating Monte Carlo estimates from this point O_i. The earliest time point is 1.
#' @param node Start generating Monte Carlo estimates from O_i node W,A or Y.
#' @param t Outcome time point of interest. 
#' @param Anode Intervention node.
#' @param intervention Specify g^*, of P(A|past).   
#' @param MC How many Monte Carlo samples should be generated.
#'
#' @return An object of class \code{tstmle}.
#' \describe{
#' \item{estimate}{Mean of the outcome at time t under specified intervention, or no intervention.}
#' \item{outcome}{Outcome at time t for each MC interation.}
#'
#' @export
#'

mcEst <- function(fit, start, node, t, Anode, intervention, MC) {

  data<-fit$data
  
  #How many in a batch:
  step<-length(grep('_0', row.names(data), value=TRUE))
  
  #Starts from time 1:
  data_lag<-fit$lag_data
  data_lag<-data.frame(data_lag[,-1])
  data_lag<-data_lag[1:(t*step),]
  estNames<-row.names(data_lag)
  
  #get actual index of Anode (based on data_lag)
  Anode<-Anode*step-1
  
  #Later on will need to include more Ws
  #This should be parallelized 
  outcome<-matrix(nrow = MC, ncol = 1)
  data_lag_origin<-data_lag
  
  #Remember the outcome for specific node being either 1 or 0
  #Need this for the clever covariate.
  res_1<-matrix(nrow=MC,ncol=1)
  res_0<-matrix(nrow=MC,ncol=1)
  
  for(B in 1:MC){
    
    data_lag<-data_lag_origin
    
    for(i in seq(start,nrow(data_lag),step)){
      
      newW<-rbinom(1,1,plogis(predict(fit$W, data_lag[i,], type="response")))
      #Update for A
      data_lag[(i+1),1]<-newW
      if(i != 1){
        #Update for A
        data_lag[(i+1),2]<-newY
      }
      
      if(is.null(intervention)){
        newA<-rbinom(1,1,plogis(predict(fit$A, data_lag[i+1,], type="response")))
        #Update for Y
        data_lag[(i+2),2]<-newW
        data_lag[(i+2),1]<-newA
      }else if(!is.null(intervention)){
        if((i+1) == Anode){
          newA<-rbinom(1,1,intervention)
          #Update for Y
          data_lag[(i+2),2]<-newW
          data_lag[(i+2),1]<-newA
        }else{
          newA<-rbinom(1,1,plogis(predict(fit$A, data_lag[i+1,], type="response")))
          #Update for Y
          data_lag[(i+2),2]<-newW
          data_lag[(i+2),1]<-newA
        }
         }

      newY<-rbinom(1,1,plogis(predict(fit$Y, data_lag[i+2,], type="response")))
      #Update for next W
      #In the last iteration this will be NA... This is ok for now, just one more row.
      data_lag[(i+3),2]<-newA
      data_lag[(i+3),1]<-newY
      
    }
    
    outcome[B,]<-newY
    
    #Get back outcome when specified node at s was either 1 or 0:
    if(node=="W"){
      if(data_lag[(start+1),1]==0){
        res_0[B,]<-newY
      }else{
        res_1[B,]<-newY
      }
    }else if(node=="A"){
      if(data_lag[(start+2),1]==0){
        res_0[B,]<-newY
      }else{
        res_1[B,]<-newY
      }
    }else if(node=="Y"){
      if(data_lag[(start+3),1]==0){
        res_0[B,]<-newY
      }else{
        res_1[B,]<-newY
      }
    }
    
  }
  
  return(list(estimate=mean(outcome), outcome=outcome, intervention=intervention,MC=MC, t=t, Anode=Anode,
              s1<-mean(res_1, na.rm=TRUE), s0=mean(res_0, na.rm = TRUE), s1_full=res_1, s0_full=res_0))

 
}

