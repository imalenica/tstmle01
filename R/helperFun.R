#####################
#Helper functions
#####################

#Get actual and predicted values for time t, based on the previously obtained fit and actual Cs.
#Minimum t is 1. 
getPred<-function(fit, t){
  
  data<-fit$data
  data_lag<-fit$lag_data
  data_lag<-data_lag[,-1]
  
  step<-length(grep('_0', row.names(data), value=TRUE))
  
  #Get actual values
  W<-data[((t*step)+1),1]
  A<-data[((t*step)+2),1]
  Y<-data[((t*step)+3),1]
  
  #Get predicted probabilities:
  W_pred<-predict(fit$W, data_lag[(t*step)-2,], type="response")
  A_pred<-predict(fit$A, data_lag[(t*step)-1,], type="response")
  Y_pred<-predict(fit$A, data_lag[(t*step),], type="response")
  
  return(list(W=W,A=A,Y=Y,W_pred=W_pred,A_pred=A_pred,Y_pred=Y_pred))
}

#Function to recalculate the EIC based on the already calculated clever covariates.
getEIC<-function(clevCov, pred_star, n){
  
  #Get all the clever covariates:
  Hy<-data.frame(clevCov$Hy)
  Ha<-data.frame(clevCov$Ha)
  Hw<-data.frame(clevCov$Hw)
  
  D<-matrix(nrow=n,ncol=1)
  
  #Calculate the EIC:
  for(i in 1:n){
    preds<-getPred(fit,i)
    D[i,]<-Hy[i,]*(preds$Y-pred_star[i,1])+Ha[i,]*(preds$A-pred_star[i,2])+Hw[i,]*(preds$W-pred_star[i,3])
  }
  
  return(Dbar=D)
}

#Sample from discrete uniform distribution
rdunif<-function(n,k) sample(1:k,n,replace=T)


