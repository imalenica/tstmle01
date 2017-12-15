#####################
# Helper functions
#####################

#' Get Predictions
#'
#' Get actual and predicted values for time t, based on the previously obtained
#' fit and actual Cs. Minimum t is 1.
#'
#' @param fit ...
#' @param t ...
#'
#' @importFrom stats predict
#
getPred <- function(fit, t) {
  data <- fit$data
  data_lag <- fit$lag_data
  data_lag <- data_lag[, -1]

  step <- length(grep("_0", row.names(data), value = TRUE))

  # Get actual values
  W <- data[((t * step) + 1), 1]
  A <- data[((t * step) + 2), 1]
  Y <- data[((t * step) + 3), 1]

  # Get predicted probabilities:
  W_pred <- stats::predict(fit$W, data_lag[(t * step) - 2, ], type = "response")
  A_pred <- stats::predict(fit$A, data_lag[(t * step) - 1, ], type = "response")
  Y_pred <- stats::predict(fit$A, data_lag[(t * step), ], type = "response")

  return(list(W = W, A = A, Y = Y, W_pred = W_pred, A_pred = A_pred,
              Y_pred = Y_pred))
}

################################################################################

#' Recalculate the EIC
#'
#' Recompute the EIF based on the previously calculated clever covariates.
#'
#' @param clevCov ...
#' @param pred_star ...
#' @param n ...
#
getEIC <- function(clevCov, pred_star, n) {

  # Get all the clever covariates:
  Hy <- data.frame(clevCov$Hy)
  Ha <- data.frame(clevCov$Ha)
  Hw <- data.frame(clevCov$Hw)

  D <- matrix(nrow = n, ncol = 1)

  # Calculate the EIC:
  for (i in seq_len(n)) {
    preds <- getPred(fit = clevCov, i)
    D[i, ] <- Hy[i, ] * (preds$Y - pred_star[i, 1]) + Ha[i, ] *
      (preds$A - pred_star[i, 2]) + Hw[i, ] * (preds$W - pred_star[i, 3])
  }
  return(Dbar = D)
}

################################################################################

#' Discrete Uniform Distribution
#'
#' Sample from discrete uniform distribution
#'
#' @param n ...
#' @param k ...
#
rdunif <- function(n, k) {
  out <- sample(seq_len(k), n, replace = TRUE)
  return(out)
}

