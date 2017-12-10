#' Main TMLE Calculations
#'
#' This function runs the TMLE and computes the estimate.
#'
#' @param clevCov Results object obtained by \code{cleverCov}. It contains clever covariates for each of 
#' the components of the likelihood, and the efficient influence curve.
#' @param fit \code{fit} object obtained by \code{initEst}.  
#'
#' @return An object of class \code{tstmle}.
#' \describe{
#' \item{fit}{Fit objects for each part of the likelihood.}
#'            
#' 
#'
#' @export
#'

mainTMLE <- function(clevCov, fit) {
  
  #Get all the clever covariates:
  Hy<-clevCov$Hy
  Ha<-clevCov$Ha
  Hw<-clevCov$Hw
  
  D<-clevCov$Dbar
  
  t<-length(Hy)
  
  Y<-matrix(nrow = t, ncol=1)
  A<-matrix(nrow = t, ncol=1)
  W<-matrix(nrow = t, ncol=1)
  
  Y_pred<-matrix(nrow = t, ncol=1)
  A_pred<-matrix(nrow = t, ncol=1)
  W_pred<-matrix(nrow = t, ncol=1)
  
  #Get actual values and predictions:
  for(i in 1:t){
    preds<-getPred(fit,i)
    
    Y[i,]<-preds$Y
    A[i,]<-preds$A
    W[i,]<-preds$W
    
    Y_pred[i,]<-preds$Y_pred
    A_pred[i,]<-preds$A_pred
    W_pred[i,]<-preds$W_pred
    
  }
  
  #TO DO: Finish until convergence
  
  if(mean(D)>sd(D)/t){
    
    #logit update:
    logUpdate_y<-glm(Y~ -1 + offset(qlogis(Y_pred)) + Hy, family='binomial')
    logUpdate_a<-glm(A~ -1 + offset(qlogis(A_pred)) + Ha, family='binomial')
    logUpdate_w<-glm(W~ -1 + offset(qlogis(W_pred)) + Hw, family='binomial')
    
    #Get epsilons
    eps_y<-logUpdate_y$coefficients
    eps_a<-logUpdate_a$coefficients
    eps_w<-logUpdate_w$coefficients
    
    #New estimate
    Y_star<-plogis(qlogis(Y_pred) + eps_y*Hy)
    A_star<-plogis(qlogis(A_pred) + eps_a*Ha)
    W_star<-plogis(qlogis(W_pred) + eps_w*Hw)
    
    Y_pred<-Y_star
    A_pred<-A_star
    W_pred<-W_star
    
  }
  
}
