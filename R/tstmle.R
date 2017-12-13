#' Time-Series Targeted Minimum Loss Estimation (tstmle)
#'
#' This function returns the final estimate of the expected value the outcome at time t, under intervention 
#' (or no intervention) as specified by the user. In addition, it returns the variance of the estimate 
#' and confidence intervals. 
#'
#' @param data data.frame object containg the time series with relevant time ordering. 
#' @param freqW A numeric specifying the Markov order for W nodes. 
#' @param freqA A numeric specifying the Markov order for A nodes. 
#' @param freqY A numeric specifying the Markov order for Y nodes. 
#' @param t Outcome time point of interest. It must be greater than the intervention node A.
#' @param Anode Intervention node.
#' @param intervention Specify g^*, of P(A|past). Right now supports only 1/0 type interventions.
#' @param B How many samples to draw from P, as part of the h-density estimation.
#' @param N How many sample to draw from P^*, as part of the h-density estimation. 
#' @param MC How many Monte Carlo samples should be generated.
#' @param maxIter Maximum number of iterations.
#' @param tol Lower bound for epsilon. 
#' @param alpha alpha
#' 
#' @return An object of class \code{tstmle}.
#' \describe{
#' \item{psi}{Estimate of the target parameter.}
#' \item{var.psi}{Variance, based on the influence function.}
#' \item{CI}{Confidence intervals.}}
#'            
#' @importFrom prodlim row.match
#'
#' @export
#'

tstmle <- function(data,freqY,freqA,freqW,t,Anode,intervention,MC=100,B=100,N=100,
                   maxIter=50,tol=10^-5,alpha=0.05) {
 
  #TO DO: add some checks for data, make sure it is in right order, format, etc.
  #TO DO: add option for different number of MCs for different parts of the estimation process.
  
  fit<-initEst(data,freqW=freqW,freqA=freqA,freqY=freqY)

  est<-mainTMLE(fit,t=t,Anode=Anode,intervention=intervention,
                B=B,N=N,MC=MC,maxIter=maxIter,tol=tol)
  
  #TO DO: add inference 
  
  
  return(list(psi=psi, var.psi=var(IC)/n), CI=list(CI_lower=CI_lower,CI_upper=CI_upper))
 
}
