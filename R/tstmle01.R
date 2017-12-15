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
#' @param intervention1 Specify g^*, of P(A|past). Right now supports only 1/0 type interventions.
#' @param intervention2 Specify g^* to compate to. Right now supports only 1/0 type interventions.
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
#' \item{IC}{Influence curve.}
#'            
#' @importFrom prodlim row.match
#'
#' @export
#'

tstmle01 <- function(data,freqY,freqA,freqW,t,Anode,intervention1,intervention2=NULL,MC=100,B=100,N=100,
                   maxIter=50,tol=10^-5,alpha=0.05) {
 
  #TO DO: add some checks for data, make sure it is in right order, format, etc.
  #TO DO: add option for different number of MCs for different parts of the estimation process.
  
  fit<-initEst(data,freqW=freqW,freqA=freqA,freqY=freqY)

  est1<-mainTMLE(fit,t=t,Anode=Anode,intervention=intervention1,
                B=B,N=N,MC=MC,maxIter=maxIter,tol=tol)
  
  est2<-mainTMLE(fit,t=t,Anode=Anode,intervention=intervention2,
                 B=B,N=N,MC=MC,maxIter=maxIter,tol=tol)
  
  diffIC<-est1$IC-est2$IC
  diffEst<-est1$psi-est2$psi
  
  #Finite sample variance of the estimator.
  var_tmle<-var(diffIC, na.rm=TRUE)/length(diffIC)
  
  #Standard error
  se_tmle<-sqrt(var_tmle)
  
  #Add CI for the difference
  ci_low<-diffEst-(stats::qnorm(1-(alpha/2)))*se_tmle
  ci_high<-diffEst+(stats::qnorm(1-(alpha/2)))*se_tmle
  
  p <-2*stats::pnorm(abs(diffEst/se_tmle), lower.tail=F)
  
  return(list(psi=diffEst, var.psi=var_tmle, CI=list(CI_lower=ci_low,CI_upper=ci_high), IC=diffIC))
 
}
