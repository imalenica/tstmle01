#' h-density estimation
#'
#' This function evaluates the h-density, necessary for the clever covariate calculation. 
#'
#' @param fit \code{fit} object obtained by \code{initEst}. 
#' @param i Where we are in the i loop (part of the clever covariate calculation).
#' @param B Number of observations to sample from P and P^*.
#' @param t Outcome time point of interest. It must be greater than the intervention node A.
#'
#' @return An object of class \code{tstmle}.
#' \describe{
#' \item{h_cy}{Empirical estimate of h_cy}
#' \item{h_ca}{Empirical estimate of h_ca}
#' \item{h_cw}{Empirical estimate of h_cw}}
#'            
#' 
#'
#' @export
#'

hEst <- function(fit, i, B, t) {
  
  #Sample B observations from P:
  p_mc<-mcEst(fit, start=1, node="A", t=t, Anode=1, MC=B, returnMC_full=TRUE)
  p<-p_mc$MCdata
  
  #How many in a batch:
  step<-length(grep('_0', row.names(p_star), value=TRUE))
  
  
  
  
  
}
