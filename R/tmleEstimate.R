#' Main TMLE Calculations
#'
#' This function runs the TMLE and computes the estimate.
#'
#' @param fit \code{fit} object obtained by \code{initEst}. 
#' @param t Outcome time point of interest. It must be greater than the intervention node A.
#' @param Anode Intervention node.
#' @param MC How many Monte Carlo samples should be generated.
#' @param maxIter Maximum number of iterations.
#' @param tol Lower bound for epsilon. 
#' 
#' @return An object of class \code{tstmle}.
#' \describe{
#' \item{psi}{Estimate of the target parameter.}
#' \item{var.psi}{Variance, based on the influence function.}}
#'            
#' 
#'
#' @export
#'

mainTMLE <- function(fit, t, Anode, MC, maxIter=50, tol=10^-5) {
  
  #Implement with one epsilon:
  iter<-0
  eps<-Inf
  
  while(iter<=maxIter & (abs(eps) > tol)){
    iter<-iter+1
    
    #Generate a new fit, using updated data:
    fit<-initEst(data, freqW = 2,freqA = 2,freqY = 2)
    clevCov<-cleverCov(fit,t=t,Anode=Anode,MC=MC)
    
    #Get all the clever covariates:
    Hy<-clevCov$Hy
    Ha<-clevCov$Ha
    Hw<-clevCov$Hw
    
    #Get the EIC:
    D<-clevCov$Dbar
    
    #Get n
    n<-length(Hy)
    
    Y<-matrix(nrow = n, ncol=1)
    A<-matrix(nrow = n, ncol=1)
    W<-matrix(nrow = n, ncol=1)
    
    Y_pred<-matrix(nrow = n, ncol=1)
    A_pred<-matrix(nrow = n, ncol=1)
    W_pred<-matrix(nrow = n, ncol=1)
    
    #Get actual values and predictions:
    for(i in 1:n){
      preds<-getPred(fit,i)
      
      Y[i,]<-preds$Y
      A[i,]<-preds$A
      W[i,]<-preds$W
      
      Y_pred[i,]<-preds$Y_pred
      A_pred[i,]<-preds$A_pred
      W_pred[i,]<-preds$W_pred
      
    }
    
    #Put all in one dataframe:
    observed<-data.frame(X=c(Y,A,W))
    pred<-data.frame(off=c(Y_pred,A_pred,W_pred))
    H<-data.frame(H=c(Hy,Ha,Hw))
    
    d<-cbind.data.frame(observed,pred,H)
    
    eps<-coef(glm(X ~ -1 + offset(qlogis(off)) + H, data=d, family='binomial'))
    eps[is.na(eps)] <- 0
    
    #Update:
    Y_star<-plogis(qlogis(Y_pred) + eps*Hy)
    A_star<-plogis(qlogis(A_pred) + eps*Ha)
    W_star<-plogis(qlogis(W_pred) + eps*Hw)
    
    #Make it easier to recombine in W,A,Y fashion.
    pred_star<-data.frame(data=c(Y_star=Y_star, A_star=A_star, W_star=W_star))
    name<-c(paste0(seq(1:n),"_C"),paste0(seq(1:n),"_B"),paste0(seq(1:n),"_A"))
    row.names(pred_star)<-name
    
    data <- data.frame(data=pred_star[order(row.names(pred_star)), ])
    
    #Generate convinient dataframe for getEIC in case this is the last iteration:
    pred_star_fin<-cbind.data.frame(Y_star=Y_star, A_star=A_star, W_star=W_star)
  
  }

  #Update EIC:
  IC<-getEIC(clevCov, pred_star_fin, n)

  return(list(psi=psi, var.psi=var(IC)/n))
    
  }
  

