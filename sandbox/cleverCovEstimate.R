#' Clever Covariate calculation
#'
#' This function calculates the clever covariate for each time point. 
#'
#' @param fit \code{fit} object obtained by \code{initEst}.  
#' @param t Outcome time point of interest. It must be greater than the intervention node A.
#' @param Anode Intervention node.
#' @param intervention Specify g^*, of P(A|past). Right now supports only 1/0 type interventions. 
#' @param MC How many Monte Carlo samples should be generated.
#'
#' @return An object of class \code{tstmle}.
#' \describe{
#' \item{fit}{Fit objects for each part of the likelihood.}
#'            
#' 
#'
#' @export
#'

cleverCov <- function(fit, t, Anode, intervention=1,  MC=1000) {

  step<-length(grep('_0', row.names(fit$data), value=TRUE))
  n<-nrow(fit$data)/step-1
  #Generate clever covariates for each of the likelihood components: W,A,Y and EIC
  
  Hy_cc<-matrix(nrow=t-1,ncol=1)
  Ha_cc<-matrix(nrow=t,ncol=1)
  Hw_cc<-matrix(nrow=t,ncol=1)
  D<-matrix(nrow=t,ncol=1)
  
  #TO DO: Probably need some kind of an internal seed for these computations.
  #Generate our P^*, intervening only on Anode. 
  p_star<-mcEst(fit, start=Anode, node="A", t=t, Anode=Anode, intervention=intervention, MC=1, returnMC=TRUE)
  
  #Add intervened data to fit:
  fit[["p_star"]]<-p_star$MCdata

  for(i in 1:n){
    
    Hy_diff<-0
    Ha_diff<-0
    Hw_diff<-0
    
    #Need to address s=t case. Think about it.
    for(s in 1:(t-1)){

      #Example for Hy:
      #Compute H_{y(s)}(C_y(i))=E[Y_{g^*}|Y(s)=1,C_y(s)=C_y(i)]-E[Y_{g^*}|Y(s)=0,C_y(s)=C_y(i)]

      if(s<i){
        #Look into the future to condition on.
        if(s==t){
          #If s=t, then just return Y^* for Hy.
          Hy_diff_add<-p_star$MCdata[length(p_star$MCdata),]
          
          Ha1<-mcEst(fit, start=s, node="Y", t=t, Anode=Anode, lag=(-1*(i-s)), MC=MC, clevCov=TRUE, set=1)
          Ha0<-mcEst(fit, start=s, node="Y", t=t, Anode=Anode, lag=(-1*(i-s)), MC=MC, clevCov=TRUE, set=0)
          Ha_diff_add<-Ha1$s-Ha0$s
          
          Hw1<-mcEst(fit, start=s, node="A", t=t, Anode=Anode, lag=(-1*(i-s)), MC=MC, clevCov=TRUE, set=1)
          Hw0<-mcEst(fit, start=s, node="A", t=t, Anode=Anode, lag=(-1*(i-s)), MC=MC, clevCov=TRUE, set=0)
          Hw_diff_add<-Hw1$s-Hw0$s 
        }else{
          Hy1<-mcEst(fit, start=s+1, node="W", t=t, Anode=Anode, lag=(-1*(i-s)), MC=MC, clevCov=TRUE, set=1)
          Hy0<-mcEst(fit, start=s+1, node="W", t=t, Anode=Anode, lag=(-1*(i-s)), MC=MC, clevCov=TRUE, set=0)
          Hy_diff_add<-Hy1$s-Hy0$s   
          
          Ha1<-mcEst(fit, start=s+1, node="Y", t=t, Anode=Anode, lag=(-1*(i-s)), MC=MC, clevCov=TRUE, set=1)
          Ha0<-mcEst(fit, start=s+1, node="Y", t=t, Anode=Anode, lag=(-1*(i-s)), MC=MC, clevCov=TRUE, set=0)
          Ha_diff_add<-Ha1$s-Ha0$s
          
          Hw1<-mcEst(fit, start=s+1, node="A", t=t, Anode=Anode, lag=(-1*(i-s)), MC=MC, clevCov=TRUE, set=1)
          Hw0<-mcEst(fit, start=s+1, node="A", t=t, Anode=Anode, lag=(-1*(i-s)), MC=MC, clevCov=TRUE, set=0)
          Hw_diff_add<-Hw1$s-Hw0$s 
        }
      }else if(s==i){
        #Usual setting.
        if(s==t){
          #If s=t, then just return Y^* for Hy.
          Hy_diff_add<-p_star$MCdata[length(p_star$MCdata),]
          
          Ha1<-mcEst(fit, start=s, node="Y", t=t, Anode=Anode, lag=0, MC=MC, clevCov=TRUE, set=1)
          Ha0<-mcEst(fit, start=s, node="Y", t=t, Anode=Anode, lag=0, MC=MC, clevCov=TRUE, set=0)
          Ha_diff_add<-Ha1$s-Ha0$s
          
          Hw1<-mcEst(fit, start=s, node="A", t=t, Anode=Anode, lag=0, MC=MC, clevCov=TRUE, set=1)
          Hw0<-mcEst(fit, start=s, node="A", t=t, Anode=Anode, lag=0, MC=MC, clevCov=TRUE, set=0)
          Hw_diff_add<-Hw1$s-Hw0$s
        }else{
          Hy1<-mcEst(fit, start=s+1, node="W", t=t, Anode=Anode, lag=0, MC=MC, clevCov=TRUE, set=1)
          Hy0<-mcEst(fit, start=s+1, node="W", t=t, Anode=Anode, lag=0, MC=MC, clevCov=TRUE, set=0)
          Hy_diff_add<-Hy1$s-Hy0$s   
          
          Ha1<-mcEst(fit, start=s+1, node="Y", t=t, Anode=Anode, lag=0, MC=MC, clevCov=TRUE, set=1)
          Ha0<-mcEst(fit, start=s+1, node="Y", t=t, Anode=Anode, lag=0, MC=MC, clevCov=TRUE, set=0)
          Ha_diff_add<-Ha1$s-Ha0$s
          
          Hw1<-mcEst(fit, start=s+1, node="A", t=t, Anode=Anode, lag=0, MC=MC, clevCov=TRUE, set=1)
          Hw0<-mcEst(fit, start=s+1, node="A", t=t, Anode=Anode, lag=0, MC=MC, clevCov=TRUE, set=0)
          Hw_diff_add<-Hw1$s-Hw0$s
        }
      }else if(s>i){
        #Look further in the past since s>i.
        if(s==t){
          #If s=t, then just return Y^* for Hy.
          Hy_diff_add<-p_star$MCdata[length(p_star$MCdata),]
          
          Ha1<-mcEst(fit, start=s, node="Y", t=t, Anode=Anode, lag=(s-i), MC=MC, clevCov=TRUE, set=1)
          Ha0<-mcEst(fit, start=s, node="Y", t=t, Anode=Anode, lag=(s-i), MC=MC, clevCov=TRUE, set=0)
          Ha_diff_add<-Ha1$s-Ha0$s 
          
          Hw1<-mcEst(fit, start=s, node="A", t=t, Anode=Anode, lag=(s-i), MC=MC, clevCov=TRUE, set=1)
          Hw0<-mcEst(fit, start=s, node="A", t=t, Anode=Anode, lag=(s-i), MC=MC, clevCov=TRUE, set=0)
          Hw_diff_add<-Hw1$s-Hw0$s
        }else{
          Hy1<-mcEst(fit, start=s+1, node="W", t=t, Anode=Anode, lag=(s-i), MC=MC, clevCov=TRUE, set=1)
          Hy0<-mcEst(fit, start=s+1, node="W", t=t, Anode=Anode, lag=(s-i), MC=MC, clevCov=TRUE, set=0)
          Hy_diff_add<-Hy1$s-Hy0$s   
          
          Ha1<-mcEst(fit, start=s+1, node="Y", t=t, Anode=Anode, lag=(s-i), MC=MC, clevCov=TRUE, set=1)
          Ha0<-mcEst(fit, start=s+1, node="Y", t=t, Anode=Anode, lag=(s-i), MC=MC, clevCov=TRUE, set=0)
          Ha_diff_add<-Ha1$s-Ha0$s 
          
          Hw1<-mcEst(fit, start=s+1, node="A", t=t, Anode=Anode, lag=(s-i), MC=MC, clevCov=TRUE, set=1)
          Hw0<-mcEst(fit, start=s+1, node="A", t=t, Anode=Anode, lag=(s-i), MC=MC, clevCov=TRUE, set=0)
          Hw_diff_add<-Hw1$s-Hw0$s
        }
      }
      
      Hy_diff<-Hy_diff+Hy_diff_add
      Ha_diff<-Ha_diff+Ha_diff_add
      Hw_diff<-Hw_diff+Hw_diff_add
    }
    
    #Store clever covariates for each time point
    Hy_cc[i,]<-Hy_diff
    Ha_cc[i,]<-Ha_diff
    Hw_cc[i,]<-Hw_diff
    
    #Calculate the EIC:
    preds<-getPred(fit,i)
    D[i,]<-Hy_cc[i,]*(preds$Y-preds$Y_pred)+Ha_cc[i,]*(preds$A-preds$A_pred)+Hw_cc[i,]*(preds$W-preds$W_pred)
    
  }
   
   return(list(Hy=Hy_cc,Ha=Ha_cc,Hw=Hw_cc, Dbar=D))
  
}