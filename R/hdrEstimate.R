#' h-density estimation
#'
#' This function evaluates the h-density ratios, necessary for the clever covariate calculation. 
#' It approximates the ratio of densities under intervened and non-intervened distributions.
#'
#' @param fit \code{fit} object obtained by \code{initEst}. 
#' @param s Where we are in the s loop (part of the clever covariate calculation).
#' @param i Where we are in the i loop (part of the clever covariate calculation).
#' @param B Number of observations to sample from P and P^*.
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

hdrEst <- function(fit, s, i, B, t, Anode, intervention=1,  MC=1000) {
  
  #Sample B observations from P:
  p_mc<-mcEst(fit, start=Anode, node="A", t=t, Anode=1, MC=B, returnMC=TRUE)
  p<-p_mc$MCdata
  
  #Sample B observations from P^*:
  p_star_mc<-mcEst(fit, start=Anode, node="A", t=t, Anode=Anode, intervention=intervention, MC=B, returnMC=TRUE)
  p_star<-p_star_mc$MCdata
  
  #Determine C based on s:

  
  
  

  }
