#' h-density estimation
#'
#' This function evaluates the h-density, necessary for the clever covariate calculation. 
#'
#' @param fit \code{fit} object obtained by \code{initEst}. 
#' @param i Where we are in the i loop (part of the clever covariate calculation).
#' @param B Number of observations to sample from P and P^*.
#' @param t Outcome time point of interest. It must be greater than the intervention node A.
#' @param Anode Intervention node.
#'
#' @return An object of class \code{tstmle}.
#' \describe{
#' \item{fit}{Fit objects for each part of the likelihood.}
#'            
#' 
#'
#' @export
#'

hEst <- function(fit, i, B, t, Anode) {
  
  #Sample B observations from P:
  p_mc<-mcEst(fit, start=Anode, node="A", t=t, Anode=1, MC=B, returnMC_full=TRUE)
  p<-p_mc$MCdata
  
  #Determine C based on s:
  
  
  
  
  
}
