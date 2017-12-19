#' Initial Estimate
#'
#' This function generates coefficients of previous variables based on a single
#' time series observation. Later on, this should include more complex
#' relationships instead of just past variables and better estimation methods.
#'
#' @param data data.frame object containg the time series with relevant time ordering.
#' @param t Outcome time point of interest. It must be greater than the intervention node A.
#' @param freqW A numeric specifying the Markov order for W nodes.
#' @param freqA A numeric specifying the Markov order for A nodes.
#' @param freqY A numeric specifying the Markov order for Y nodes.
#'
#' @return An object of class \code{tstmle01}.
#' \describe{
#' \item{fitW}{Fit object for W part of the likelihood.}
#' \item{fitA}{Fit object for A part of the likelihood.}
#' \item{fitY}{Fit object for Y part of the likelihood.}
#' \item{probW}{Probability of success for the W process.} 
#' \item{probA}{Probability of success for the A process.}
#' \item{probY}{Probability of success for the Y process.}  
#' \item{freqW}{Input Markov order for W.}
#' \item{freqA}{Input Markov order for A.}
#' \item{freqY}{Input Markov order for Y.}
#' \item{data}{Input data.}
#' \item{lag_data}{Data with lagged components.}
#' }
#'
#' @importFrom Hmisc Lag
#' @importFrom stats glm as.formula quasibinomial complete.cases
#'
#' @export
#

initEst <- function(data, t, freqW = NULL, freqA = NULL, freqY = NULL) {
  
  if (all.equal(freqA, freqW, freqY)) {

    # How many in a batch:
    step <- length(grep("_1$", row.names(data), value = TRUE))
    
    #Note: time-series needs to be of length n+t for the estimation to make sense for later lags. 
    #Considered time-series might need to be longer that what we need to get it. 
    #Ex: for size n=100 and t=5, we need a time-series of 105 points. 
    n <- (nrow(data) / step - 1) - (t-1)
    n_true <- (nrow(data) / step - 1)
    
    data_orig<-data
    
    names<-row.names(data)
    data<-data.frame(data=data[1:((step*n)+step),])
    row.names(data)<-names[1:((step*n)+step)]
    
    # Lag past
    res <- lapply(seq_len(freqW), function(x) {
      Lag(data[, 1], x)
    })
    res <- data.frame(matrix(unlist(res), nrow = length(res[[1]])),
                      stringsAsFactors = FALSE)

    data_est_full <- cbind.data.frame(data = data, res)

    # Drop time 0 for estimation
    cc <- stats::complete.cases(data_est_full)
    data_est <- data_est_full[cc, ]
    data_est <- data.frame(data_est[-1, ])

    # Separate by process:
    W <- data_est[grep("W", row.names(data_est), value = TRUE), ]
    A <- data_est[grep("A", row.names(data_est), value = TRUE), ]
    Y <- data_est[grep("Y", row.names(data_est), value = TRUE), ]

    # Estimate the process coefficients. Should allow better strategies later on
    fitW <- stats::glm(formula = stats::as.formula("data ~ ."),
                       family = stats::quasibinomial(link = "logit"), data = W)
    fitA <- stats::glm(formula = stats::as.formula("data ~ ."),
                       family = stats::quasibinomial(link = "logit"), data = A)
    fitY <- stats::glm(formula = stats::as.formula("data ~ ."),
                       family = stats::quasibinomial(link = "logit"), data = Y)
    
    W_est<-data.frame(W[,-1])
    A_est<-data.frame(A[,-1])
    Y_est<-data.frame(Y[,-1])
    
    #Get probabilities of success:
    probW<-data.frame(prob=stats::predict(fitW,W_est,type = "response"))
    probA<-data.frame(prob=stats::predict(fitA,A_est,type = "response"))  
    probY<-data.frame(prob=stats::predict(fitY,Y_est,type = "response"))  
    
    #TO DO: Cases where we don't have enough data to estimate all probabilities
    #Create a file that has probabilities and lag combinations:
    W_comb<-cbind.data.frame(W_est,probW)
    A_comb<-cbind.data.frame(A_est,probA)
    Y_comb<-cbind.data.frame(Y_est,probY)
      
    unique_W_comb<-unique(W_comb)
    row.names(unique_W_comb)<-NULL
    unique_A_comb<-unique(A_comb)
    row.names(unique_A_comb)<-NULL
    unique_Y_comb<-unique(Y_comb)
    row.names(unique_Y_comb)<-NULL
    
    fit <- list(W = fitW, A = fitA, Y = fitY, 
                estW=W_est, estA=A_est, estY=Y_est,
                combW=unique_W_comb,combA=unique_A_comb,combY=unique_Y_comb,
                freqW = freqW, freqA = freqA, freqY = freqY,
                n=n, n_true=n_true, step=step,
                data = data, data_orig=data_orig, 
                lag_data = data_est)
  }
  return(fit)
}

