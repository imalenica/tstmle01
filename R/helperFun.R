#####################
#Helper functions
#####################

#Get actual and predicted values for time t, based on the previously obtained fit and actual Cs.
#Minimum t is 1. 
getPred<-function(fit, t){
  
  data<-fit$data
  data_lag<-fit$lag_data
  
  step<-length(grep('_0', row.names(data), value=TRUE))
  
  #Get actual values
  W<-data[((t*step)+1),1]
  A<-data[((t*step)+2),1]
  Y<-data[((t*step)+3),1]
  
  #Get predicted values
  W_pred<-predict(fit$W, data_lag[(t*step)-2,], type="response")
  A_pred<-predict(fit$A, data_lag[(t*step)-1,], type="response")
  Y_pred<-predict(fit$A, data_lag[(t*step),], type="response")
  
  return(list(W=W,A=A,Y=Y,W_pred=W_pred,A_pred=A_pred,Y_pred=Y_pred))
}

