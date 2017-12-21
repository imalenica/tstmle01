#' Main TMLE Calculations
#'
#' This function runs the TMLE and computes the estimate.
#'
#' @param fit \code{fit} object obtained by \code{initEst}.
#' @param t Outcome time point of interest. It must be greater than the intervention node A.
#' @param Anode Intervention node.
#' @param intervention Specify %g^*, of %P(A \mid \text{past}). Right now, this supports only 1/0 type interventions.
#' @param B How many samples to draw from P, as part of the h-density estimation.
#' @param N How many sample to draw from %P^*, as part of the h-density estimation.
#' @param MC How many Monte Carlo samples should be generated.
#' @param maxIter Maximum number of iterations.
#' @param tol Lower bound for epsilon.
#' @param alpha alpha
#'
#' @return An object of class \code{tstmle01}.
#' \describe{
#' \item{psi}{Estimate of the target parameter.}
#' \item{var.psi}{Variance, based on the influence function.}}
#' \item{IC}{Influence curve.}
#'
#' @importFrom stats var qnorm plogis qlogis glm coef
#'
#' @export
#

mainTMLE <- function(fit, t, Anode, intervention = NULL, alpha = 0.05, B = 100,
                     N = 100, MC = 100, maxIter = 50, tol=10 ^ -3) {

  # How many in a batch:
  step <- fit$step
  
  # Get n
  n <- fit$n
  
  #Actual n
  n_true <- fit$n_true
  
  # TO DO: Probably need some kind of an internal seed for these computations.
  # Generate our P^*, intervening only on Anode, from Anode.
  p_star <- mcEst(fit, start = Anode, node = "A", t = t, Anode = Anode,
                  intervention = intervention, MC = 1, returnMC_full = TRUE, full=TRUE)
  fit[["p_star"]] <- p_star$MCdata
  
  # Sample N observations from P^*:
  # Need to sample the full time-series because of the i-th comparison.
  p_star_mc <- mcEst(fit, start = 1, node = "W", t = t, Anode = Anode,
                     intervention = intervention, MC = N, returnMC_full = TRUE, full=TRUE)
  fit[["h_star"]] <- p_star_mc$MCdata
  
  # Sample B observations from P:
  p_mc <- mcEst(fit, start = 1, node = "A", t = t, Anode = 1, MC = B, returnMC_full = TRUE, full=TRUE)
  fit[["h"]] <- p_mc$MCdata
  
  ###################################
  # Update
  ##################################
  
  # Implement with one epsilon:
  iter <- 0
  eps <- Inf
  
  while (iter <= maxIter & (abs(eps) > tol)) {
 
    iter <- iter + 1
    
    if(iter==1){
      # Calculate the clever covariate with the initial fit
      clevCov <- cleverCov(fit, t = t, Anode = Anode, intervention = intervention, 
                           B = B, N = N, MC = MC)
    }else{
      # Calculate the clever covariate with the new fit
      clevCov <- cleverCov(fit, t = t, Anode = Anode, intervention = intervention, 
                           B = B, N = N, MC = MC, update=TRUE)
    }

    # Get all the clever covariates:
    Hy <- clevCov$Hy
    Ha <- clevCov$Ha
    Hw <- clevCov$Hw

    # Get the non-centered EIC:
    D <- clevCov$Dbar
    
    Y <- matrix(nrow = n, ncol = 1)
    A <- matrix(nrow = n, ncol = 1)
    W <- matrix(nrow = n, ncol = 1)

    Y_pred <- matrix(nrow = n, ncol = 1)
    A_pred <- matrix(nrow = n, ncol = 1)
    W_pred <- matrix(nrow = n, ncol = 1)

    # Get actual values and predictions:
    for (i in seq_len(n)) {
      
      preds <- getPred(fit, i)

      Y[i, ] <- preds$Y
      A[i, ] <- preds$A
      W[i, ] <- preds$W

      Y_pred[i, ] <- preds$Y_pred
      A_pred[i, ] <- preds$A_pred
      W_pred[i, ] <- preds$W_pred
    }

    # Put all in one dataframe:
    observed <- data.frame(X = c(Y, A, W))
    if(iter==1){
      pred <- data.frame(off = c(Y_pred, A_pred, W_pred))
    }else{
      #Returns log odds, essentially.
      Y_pred<-exp(Y_pred)/(1+exp(Y_pred))
      A_pred<-exp(A_pred)/(1+exp(A_pred))
      W_pred<-exp(W_pred)/(1+exp(W_pred))

      pred <- data.frame(off = c(Y_pred, A_pred, W_pred))
    }
    
    H <- data.frame(H = c(Hy, Ha, Hw))

    d <- cbind.data.frame(observed, pred, H)

    #off is prob
    eps <- stats::coef(stats::glm(X ~ -1 + stats::offset(stats::qlogis(off)) + H, data = d, family = "quasibinomial"))[2]
    eps[is.na(eps)] <- 0

    # Update (probabilities):
    Y_star <- stats::plogis(stats::qlogis(Y_pred) + eps * Hy)
    A_star <- stats::plogis(stats::qlogis(A_pred) + eps * Ha)
    W_star <- stats::plogis(stats::qlogis(W_pred) + eps * Hw)
    
    #Create new unique combinations. Not really necessary, mostly checks:
    
    #Create a file that has probabilities and lag combinations:
    W_comb<-cbind.data.frame(fit$estW,W_star)
    A_comb<-cbind.data.frame(fit$estA,A_star)
    Y_comb<-cbind.data.frame(fit$estY,Y_star)
    
    unique_W_comb<-unique(W_comb)
    row.names(unique_W_comb)<-NULL
    unique_A_comb<-unique(A_comb)
    row.names(unique_A_comb)<-NULL
    unique_Y_comb<-unique(Y_comb)
    row.names(unique_Y_comb)<-NULL

    #Update log odds
    Y_odd <- data.frame(data=stats::qlogis(Y_pred) + eps * Hy)
    A_odd <- data.frame(data=stats::qlogis(A_pred) + eps * Ha)
    W_odd <- data.frame(data=stats::qlogis(W_pred) + eps * Hw)
    
    row.names(Y_odd)<-paste0("Y_",1:nrow(Y_odd))
    row.names(A_odd)<-paste0("A_",1:nrow(A_odd))
    row.names(W_odd)<-paste0("W_",1:nrow(W_odd))
    
    #Get new coeffs:
    W_odd<-cbind.data.frame(W_odd,fit$estW)
    A_odd<-cbind.data.frame(A_odd,fit$estA)
    Y_odd<-cbind.data.frame(Y_odd,fit$estY)
    
    fit$W <- stats::glm(formula = stats::as.formula("data ~ ."), data = W_odd)
    fit$A <- stats::glm(formula = stats::as.formula("data ~ ."), data = A_odd)
    fit$Y <- stats::glm(formula = stats::as.formula("data ~ ."), data = Y_odd)
    
  }

  # Get final estimate at time t, for a parameter that does not average across batches.
  # TO DO: implement the J-type parameter. Better for actually stationary data!
  psi <- data[(step * (t + 1)), ]
  Psi<-rep(psi,nrow(data))
  
  # Update EIC:
  IC <- clevCov$Dbar - Psi

  # Finite sample variance of the estimator.
  var_tmle <- stats::var(IC, na.rm = TRUE) / length(IC)

  # Standard error
  se_tmle <- sqrt(var_tmle)

  # Add CI for estimate
  ci_low <- psi - (stats::qnorm(1 - (alpha / 2))) * se_tmle
  ci_high <- psi + (stats::qnorm(1 - (alpha / 2))) * se_tmle

  return(list(psi = psi, var.psi = var_tmle,
              CI = list(CI_lower = ci_low, CI_upper = ci_high), IC = IC),
              conv = abs(epsilon) <= tol, iter=iter)
}

