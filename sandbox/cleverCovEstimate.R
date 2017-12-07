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
  
  #Generate clever covariates for each of the likelihood components: W,A,Y
  
  Hy<-matrix(nrow=t,ncol=1)
  Ha<-matrix(nrow=t,ncol=1)
  Hw<-matrix(nrow=t,ncol=1)
  
  #TO DO: Probably need some kind of an internal seed for these computations.
  #Generate our P^*, intervening only on Anode. 
  p_star<-mcEst(fit, start=Anode, node="A", t=t, Anode=Anode, intervention=intervention, MC=1, init=TRUE, returnMC=TRUE)
  
  #Add intervened data to fit:
  fit[["p_star"]]<-p_star$MCdata
  
  for(i in 1:t){
    
    Hy_diff<-0
    
    for(s in 1:t){

      #Q: What if, for example we start generating MCs from node 4, but our intervention is at 3? 
      #A: Perhaps our input data for this step should be an MC with intervention already there?
      
      #Compute H_{y(s)}(C_y(i))=E[Y_{g^*}|Y(s)=1,C_y(s)=C_y(i)]-E[Y_{g^*}|Y(s)=0,C_y(s)=C_y(i)]
      #For H_y, start node is W.
      if(s<i){
        if(s==t){
          #If s=t, then just return E[Y^*=1|C_y]-E[Y^*=0|C_y]
          Hy<-mcEst(fit, start=, node="W", t=t, Anode=Anode, intervention=intervention, lag=(-1*(i-s)), MC=MC)
          Hy_diff_add<-Hy$s1-Hy$s0
        }else{
          Hy<-mcEst(fit, start=s+1, node="W", t=t, Anode=Anode, intervention=intervention, lag=(-1*(i-s)), MC=MC)
          Hy_diff_add<-Hy$s1-Hy$s0          
        }
      }else if(s==i){
        if(s==t){
          #If s=t, then just return E[Y^*=1|C_y]-E[Y^*=0|C_y]
          Hy<-mcEst(fit, start=, node="W", t=t, Anode=Anode, intervention=intervention, lag=(-1*(i-s)), MC=MC)
          Hy_diff_add<-Hy$s1-Hy$s0
        }else{
          Hy<-mcEst(fit, start=s+1, node="W", t=t, Anode=Anode, lag=0, MC=MC,clevCov=TRUE)
          Hy_diff_add<-Hy$s1-Hy$s0         
        }
      }else if(s>i){
        if(s==t){
          #If s=t, then just return E[Y^*=1|C_y]-E[Y^*=0|C_y]
          Hy<-mcEst(fit, start=, node="W", t=t, Anode=Anode, intervention=intervention, lag=(-1*(i-s)), MC=MC)
          Hy_diff_add<-Hy$s1-Hy$s0
        }else{
          Hy<-mcEst(fit, start=s+1, node="W", t=t, Anode=Anode, intervention=intervention, lag=(s-i), MC=MC)
          Hy_diff_add<-Hy$s1-Hy$s0       
        }
      }
      
      Hy_diff<-Hy_diff+Hy_diff_add
      
      
    }
    
    
    
  }
  
  mcEst(fit, start=1, node="W", t=5, Anode=3, intervention=1, lag=1, MC=1000)
  
  
  
  
  
  
  
}