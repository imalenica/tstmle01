#' Compute the Monte Carlo estimate
#'
#' This function computes the Monte Carlo estimate of the expected value of the parameter of interest
#' under specified intervention and under the estimated mechanism. 
#' 
#' @param fit \code{fit} object obtained by \code{initEst}. 
#' @param start Start generating Monte Carlo estimates from this point O_i. The earliest time point is 1. For the clever covariate, this is also our s.
#' @param node Start generating Monte Carlo estimates from O_i node W,A or Y.
#' @param t Outcome time point of interest. It must be greater than the intervention node A.
#' @param Anode Intervention node.
#' @param intervention Specify g^*, of P(A|past). Right now supports only 1/0 type interventions. 
#' @param lag This is an user impossed dependency lag necessary for the calculation of the clever covariate in the targeting step. 
#' It refers to O_i, and it is of the same dimension as the actual C(i), but further in the past.
#' Default is 0, meaning that it starts right after the node in question. Also, note that lag should be equal or smaller than start.   
#' @param MC How many Monte Carlo samples should be generated.
#' @param init If TRUE, estimates the mean difference between the estimated Y_g*= 1 and 
#' Y_g* = 0 (Y under intervention). If FALSE, we are within a s loop and look at the Y before starting MC draws. 
#' @param returnMC If TRUE, returns all MC draws. 
#' @param clevCov If TRUE, this MC is used for the calculation of the clever covariate. Instead of observed data,
#' it used intervened P* for further MC draws.
#'
#' @return An object of class \code{tstmle}.
#' \describe{
#' \item{estimate}{Mean of the outcome at time t under specified intervention, or no intervention.}
#' \item{outcome}{Outcome at time t for each MC interation.}
#' \item{intervention}{Intervention specified.}
#' \item{MC}{How many Monte Carlo samples should be generated.}
#' \item{Anode}{Intervention node as a function of O_i.}
#' \item{s1}{Mean of the intervened outcome given s=1 (used for the clever covariate calculation).} 
#' \item{s0}{Mean of the intervened outcome given s=0 (used for the clever covariate calculation).} 
#' \item{s1_full}{Intervened outcomes given s=1.}  
#' \item{s0_full}{Intervened outcomes given s=0.}  
#' \item{MCdata}{If \code{returnMC} is TRUE, returns a data.frame with MC time-series.} 
#' 
#' @export
#'

mcEst <- function(fit, start=1, node="W", t, Anode, intervention=NULL, lag=0, MC, init=FALSE, returnMC=FALSE, clevCov=FALSE) {

  if(clevCov==TRUE){
    data<-fit$p_star
  }else{
    data<-fit$data
  }

  #How many in a batch:
  step<-length(grep('_0', row.names(data), value=TRUE))

  #Impose artifical order (C(i)) for the clever covariate calculation. 
  #Should default to whatever was estimated in initEst()
  
  if(lag<0){
    #If lag is negative, move the time series up w.r.t. reference
    
    #This gives further up lag first
    res<-lapply((step*lag+1):(fit$freqW+step*lag), function(x) {Lag(data[,1], x)})
    res <- data.frame(matrix(unlist(res), nrow=length(res[[1]])), stringsAsFactors = FALSE)
    
    #TO DO: Make this more general later on. To time right now
    res<-res[,c(2,1)]
    
    data_est_full<-cbind.data.frame(data=data,res)
    
    #Drop time 0 for estimation
    cc<-complete.cases(data_est_full)
    data_est<-data_est_full[cc,]
    data_lag<-data.frame(data_est[-1,])
    data_lag<-data.frame(data_lag[,-1])
    
    #Now have to take into account it does not start at 1...
    if(clevCov==FALSE){
      data_lag<-data_lag[1:((t*step)-(lag*step)),]
      estNames<-row.names(data_lag)  
    }
  }else if(lag>=0){
    
    #This is ok under the assumption orders are the same dimension.
    #Note that "lags" in the clever covariate are defined by the size of the "step".
    #TO DO: Change this to an option where we can have varying dimensions for W,A,Y.
    
    res<-lapply((1+step*lag):(fit$freqW+step*lag), function(x) {Lag(data[,1], x)})
    res <- data.frame(matrix(unlist(res), nrow=length(res[[1]])), stringsAsFactors = FALSE)
    
    data_est_full<-cbind.data.frame(data=data,res)
    
    #Drop time 0 for estimation
    #Notice: drop step*(lag+1) each time
    cc<-complete.cases(data_est_full)
    data_est<-data_est_full[cc,]
    data_lag<-data.frame(data_est[-1,])
    data_lag<-data.frame(data_lag[,-1])
    
    #Now have to take into account it does not start at 1...
    if(clevCov==FALSE){
      data_lag<-data_lag[1:((t*step)-(lag*step)),]
      estNames<-row.names(data_lag)  
    }
  }
  
  #get actual index of Anode 
  Anode<-Anode*step+2-(lag+1)*step
  
  #Get actual index for start, depended on W.
  #(right now, start is set to be the batch, not actual index)
  start<-start*step-(lag+1)*step+1
  
  #TO DO: Later on will need to include more Ws
  #TO DO: This should be parallelized 
  outcome<-matrix(nrow = MC, ncol = 1)
  data_lag_origin<-data_lag
  
  #Remember the outcome for specific node being either 1 or 0
  #Need this for the clever covariate.
  res_1<-matrix(nrow=MC,ncol=1)
  res_0<-matrix(nrow=MC,ncol=1)
  
  #Prepare to return all MCs (up to time t generated draws)
  retMC<-matrix(nrow = (nrow(data_lag)+step), ncol=MC)
  row.names(retMC)<-row.names(data)[1:((t+1)*step)] 

  for(B in 1:MC){
    
    data_lag<-data_lag_origin
    
    #Keep track where we are in the loop
    iter<-1
    
    #i always dependent on W iterations. Make everything offset of W
    for(i in seq(start,nrow(data_lag),step)){
      
      #Possibly need to skip the update of W and A at the first iteration, 
      #depending on from which node we start. Hence the need for iter and ommitting some nodes.
      
      #Possibly need to condition on earlier time points for calculating our clever covariate. 
      #This is our i in the clever covariate. 
      
      if(iter==1){
        
        if(node=="A"){
          
          #Start at A.
          if(is.null(intervention)){
            newA<-rbinom(1,1,plogis(predict(fit$A, data_lag[i+1,], type="response")))
            #Update for Y
            data_lag[(i+2),1]<-newA
          }else if(!is.null(intervention)){
            if((i+1) == Anode){
              newA<-rbinom(1,1,intervention)
              #Update for Y
              data_lag[(i+2),1]<-newA
            }else{
              newA<-rbinom(1,1,plogis(predict(fit$A, data_lag[i+1,], type="response")))
              #Update for Y
              data_lag[(i+2),1]<-newA
            }
          }
          
          newY<-rbinom(1,1,plogis(predict(fit$Y, data_lag[i+2,], type="response")))
          #Update for next W
          data_lag[(i+3),2]<-newA
          data_lag[(i+3),1]<-newY
          
          iter<-iter+1
          
        }else if(node=="Y"){
          
          #Now we skip both W and A, start generating MC draws from Y.
          newY<-rbinom(1,1,plogis(predict(fit$Y, data_lag[i+2,], type="response")))
          #Update for next W
          data_lag[(i+3),1]<-newY
          
          iter<-iter+1
        }else if(node=="W"){
          
          #Don't skip anything, start with W.
          newW<-rbinom(1,1,plogis(predict(fit$W, data_lag[i,], type="response")))
          #Update for A
          data_lag[(i+1),1]<-newW
          
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
          
          iter<-iter+1
          
        }
        }else{
          
          #For later iterations, do as always.
          #MC draws start from "start" and W node.
          newW<-rbinom(1,1,plogis(predict(fit$W, data_lag[i,], type="response")))
          #Update for A
          data_lag[(i+1),1]<-newW
          data_lag[(i+1),2]<-newY
          
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
          
          iter<-iter+1
          
        }
    }
    
    outcome[B,]<-newY
    
    if(returnMC==TRUE){
      data_lag<-data_lag[-1,1]
      #Collect changed part of the series (from start until the end)
      if(node=="W"){
        retMC[((start+step):nrow(retMC)),B]<-data_lag[(start:length(data_lag))]
        retMC[(1:(start+step-1)),B]<-data[(1:(start+step-1)),1]
      }else if(node=="A"){
        retMC[((start+step+1):nrow(retMC)),B]<-data_lag[((start+1):length(data_lag))]
        retMC[(1:(start+step)),B]<-data[(1:(start+step)),1]
      }else if(node=="Y"){
        retMC[((start+step+2):nrow(retMC)),B]<-data_lag[((start+2):length(data_lag))]
        retMC[(1:(start+step+1)),B]<-data[(1:(start+step+1)),1]
      }
      
    }
    
    #Get back outcome when specified node at s (or Y^*) was either 1 or 0:
    if(init=="TRUE"){
      #If we are not in the s loop, and want to see the difference between newY. (Y^*)
      if(newY==0){
        res_0[B,]<-newY
      }else{
        res_1[B,]<-newY
      }
    }else{
      #Look at the node at s (it will be the one before start)
      if(data_lag[(start-1),1]==0){
        res_0[B,]<-newY
      }else{
        res_1[B,]<-newY
      }
    }

  }

  if(returnMC==TRUE){
    return(list(estimate=mean(outcome), outcome=outcome, intervention=intervention,MC=MC, t=t, Anode=Anode,
                s1=mean(res_1, na.rm=TRUE), s0=mean(res_0, na.rm = TRUE), s1_full=res_1, s0_full=res_0,MCdata=retMC))
    
  }else{
    return(list(estimate=mean(outcome), outcome=outcome, intervention=intervention,MC=MC, t=t, Anode=Anode,
                s1=mean(res_1, na.rm=TRUE), s0=mean(res_0, na.rm = TRUE), s1_full=res_1, s0_full=res_0))
    
  }

}

