#' h^*-density estimation
#'
#' This function evaluates the h^*-density, necessary for the clever covariate calculation. 
#'
#' @param fit \code{fit} object obtained by \code{initEst}. 
#' @param s Where we are in the s loop (part of the clever covariate calculation).
#' @param i Where we are in the i loop (part of the clever covariate calculation).
#' @param B Number of observations to sample from P and P^*.
#' @param t Outcome time point of interest. It must be greater than the intervention node A.
#' @param Anode Intervention node.
#' @param intervention Specify g^*, of P(A|past). Right now supports only 1/0 type interventions. 
#'
#' @return An object of class \code{tstmle}.
#' \describe{
#' \item{h_cy_star}{Empirical estimate of h_cy under the intervention g^*.}
#' \item{h_ca_star}{Empirical estimate of h_ca under the intervention g^*.}
#' \item{h_cw_star}{Empirical estimate of h_cw under the intervention g^*.}}
#'            
#' 
#'
#' @export
#'

h_starEst <- function(fit, s, i, B, t, Anode, intervention=1) {

  #Sample B observations from P^*:
  #Need to sample the full time-series because of the i-th comparison. 
  p_star_mc<-mcEst(fit, start=1, node="A", t=t, Anode=Anode, intervention=intervention, MC=B, returnMC_full=TRUE)
  p_star<-p_star_mc$MCdata
  
  #How many in a batch:
  step<-length(grep('_0', row.names(p_star), value=TRUE))

  #Determine Cs based on s:
  #Once again, we assume the same Markov order for W,A,Y.
  #TO DO: Change this assumption in future implementation
  
  Cy<-data.frame(matrix(nrow=B,ncol=fit$freqW))
  Ca<-data.frame(matrix(nrow=B,ncol=fit$freqW))
  Cw<-data.frame(matrix(nrow=B,ncol=fit$freqW))
  
  Cy_i<-data.frame(matrix(nrow=B,ncol=fit$freqW))
  Ca_i<-data.frame(matrix(nrow=B,ncol=fit$freqW))
  Cw_i<-data.frame(matrix(nrow=B,ncol=fit$freqW))
  
  Cy_match<-data.frame(matrix(nrow=B,ncol=1))
  Ca_match<-data.frame(matrix(nrow=B,ncol=1))
  Cw_match<-data.frame(matrix(nrow=B,ncol=1))
  
  for(b in 1:B){
    res<-lapply(1:fit$freqW, function(x) {Lag(p_star[,b], x)})
    res <- data.frame(matrix(unlist(res), nrow=length(res[[1]])), stringsAsFactors = FALSE)
    data_est_full<-cbind.data.frame(data=p_star[,b],res)
    data_est_full<-data_est_full[,-1]

    #Which C_y(s)/C_a(s)/C_w(s) do we care about:
    Cy[b,]<-data_est_full[(step+s*step),]
    Ca[b,]<-data_est_full[(step+s*step-1),]
    Cw[b,]<-data_est_full[(step+s*step-2),]  
    
    #Which C_y(i)/C_a(i)/C_w(i) do we care about:
    Cy_i[b,]<-data_est_full[(step+i*step),]
    Ca_i[b,]<-data_est_full[(step+i*step-1),]
    Cw_i[b,]<-data_est_full[(step+i*step-2),]  
    
    Cy_match[b,]<-prodlim::row.match(Cy[b,],Cy_i[b,])
    Ca_match[b,]<-prodlim::row.match(Ca[b,],Ca_i[b,])
    Cw_match[b,]<-prodlim::row.match(Cw[b,],Cw_i[b,])
    
  }
  
  h_cy<-sum(Cy_match, na.rm = TRUE)/B
  h_ca<-sum(Ca_match, na.rm = TRUE)/B
  h_cw<-sum(Cw_match, na.rm = TRUE)/B
  
  return(list(h_cy_star=h_cy,h_ca_star=h_ca,h_cw_star=h_cw))
  }
