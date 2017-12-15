#' Compute the Monte Carlo estimate
#'
#' This function computes the Monte Carlo estimate of the expected value of the
#' parameter of interest under specified intervention and under the estimated
#' mechanism.
#'
#' @param fit \code{fit} object obtained by \code{initEst}.
#' @param start Start generating Monte Carlo estimates from this point %O_i. The
#'  earliest time point is 1. For the clever covariate, this is also our s.
#' @param node Start generating Monte Carlo estimates from O_i node W, A or Y.
#' @param t Outcome time point of interest. It must be greater than the
#'  intervention node A.
#' @param Anode Intervention node.
#' @param intervention Specify %g^*, of %P(A \mid \text{past}). Right now, this
#'  supports only 1/0 type interventions.
#' @param lag This is an user impossed dependency lag necessary for the
#'  calculation of the clever covariate in the targeting step. It refers to
#'  %O_i, and it is of the same dimension as the actual C(i), but further in the
#'  past. Default is 0, meaning that it starts right after the node in question.
#'  Also, note that lag should be equal or smaller than start.
#' @param MC How many Monte Carlo samples should be generated.
#' @param returnMC If \code{TRUE}, returns all MC draws up until outcome time.
#' @param returnMC_full If \code{TRUE}, returns all full MC draws.
#' @param clevCov If \code{TRUE}, this MC is used for the calculation of the
#'  clever covariate. Instead of observed data, it used intervened %P^* for
#'  further MC draws.
#' @param set Set the s node to either 1 or 0. Used for the clever covariate
#'  calculation.
#'
#' @return An object of class \code{tstmle}.
#' \describe{
#' \item{estimate}{Mean of the outcome at time t under specified intervention,
#'  or no intervention.}
#' \item{outcome}{Outcome at time t for each MC interation.}
#' \item{intervention}{Intervention specified.}
#' \item{MC}{How many Monte Carlo samples should be generated.}
#' \item{Anode}{Intervention node as a function of O_i.}
#' \item{s}{Mean of the intervened outcome given s=1 or s=0 (used for the clever
#'  covariate calculation).}
#' \item{s_full}{Intervened outcomes given s=1 or s=0.}
#' \item{MCdata}{If \code{returnMC} is \code{TRUE}, returns a \code{data.frame}
#'  with MC time-series.}
#' }
#'
#' @importFrom stats complete.cases rbinom plogis predict
#'
#' @export
#
mcEst <- function(fit, start = 1, node = "W", t, Anode, intervention = NULL,
                  lag = 0, MC, returnMC = FALSE, returnMC_full = FALSE,
                  clevCov = FALSE, set = NULL) {

  # Checks
  if (start < lag) {
    stop("Start time has to be equal or greater than how far in the past we want to go. Pick some start>=lag.")
  }

  if (start <= 0) {
    stop("First point in a time series is 1.")
  }

  if (Anode <= 0) {
    stop("Anode must be a positive number indicating O_i in a time series.")
  }

  if (t <= 0) {
    stop("Outcome must be a positive number indicating O_i in a time series.")
  }

  if (t < Anode) {
    warning("Still intervention will not have any effect on the specified outcome.")
  }

  # How many in a batch:
  step <- length(grep("_0", row.names(data), value = TRUE))

  # If estMC is used for clever covariate calculation, use P^*.
  if (clevCov == TRUE) {
    data <- fit$p_star

    # Set s node to "set" (1/0), for Y, A, W.
    if (node == "W") {
      s <- start
      data[(s * step), ] <- set
    } else if (node == "Y") {
      s <- start
      data[(s * step) + 2, ] <- set
    } else if (node == "A") {
      s <- start
      data[(s * step) + 1, ] <- set
    }
  } else {
    data <- fit$data
  }

  # Prepare to return all MCs (up to time t (outcome) generated draws)
  if (returnMC == returnMC_full) {
    warning("Can't return both a full time-series and a shorter time-series. Returning full time-series.")
  }

  if (returnMC == TRUE) {
    retMC <- matrix(nrow = (t * step + step), ncol = MC)
    row.names(retMC) <- row.names(data)[1:((t + 1) * step)]
  }

  if (returnMC_full == TRUE) {
    retMC <- matrix(nrow = nrow(data), ncol = MC)
    row.names(retMC) <- row.names(data)
  }

  # Impose artifical order (C(i)) for the clever covariate calculation.
  # Should default to whatever was estimated in initEst()

  if (lag < 0) {
    # If lag is negative, move the time series up w.r.t. reference

    # This gives further up lag first
    res <- lapply((step * lag + 1):(fit$freqW + step * lag), function(x) {
      Lag(data[, 1], x)
    })
    res <- data.frame(matrix(unlist(res), nrow = length(res[[1]])),
                      stringsAsFactors = FALSE)

    data_est_full <- cbind.data.frame(data = data, res)
    data_lag <- data.frame(data_est_full[, -1])

    n <- nrow(data_lag) / step - 1

    # Now have to take into account it does not start at 1...
    if (clevCov == FALSE & returnMC_full == FALSE) {
      data_lag <- data_lag[1:((t * step) - (lag * step)), ]
      estNames <- row.names(data_lag)
    }

    # Now that we have all future lags, shorten the time-series.
    if (returnMC_full == FALSE & n > t) {
      data_lag <- data_lag[1:((t + 1) * step), ]
      estNames <- row.names(data_lag)[1:((t + 1) * step)]
    }

    # get actual index of Anode
    Anode <- Anode * step + 2

    # Get actual index for start, depended on W.
    # (right now, start is set to be the batch, not actual index)
    start <- start * step + 1
  } else if (lag >= 0) {

    # This is ok under the assumption orders are the same dimension.
    # Note that "lags" in the clever covariate are defined by the size of the "step".
    # TODO: Change this to an option where we can have varying dimensions for W,A,Y.

    res <- lapply((1 + step * lag):(fit$freqW + step * lag), function(x) {
      Lag(data[, 1], x)
    })
    res <- data.frame(matrix(unlist(res), nrow = length(res[[1]])),
                      stringsAsFactors = FALSE)

    data_est_full <- cbind.data.frame(data = data, res)

    # Drop time 0 for estimation
    # Notice: drop step*(lag+1) each time
    cc <- stats::complete.cases(data_est_full)
    data_est <- data_est_full[cc, ]
    data_lag <- data.frame(data_est[-1, ])
    data_lag <- data.frame(data_lag[, -1])

    n <- nrow(data_lag) / step

    # Option for returnMC
    if (clevCov == FALSE & returnMC_full == FALSE) {
      data_lag <- data_lag[1:((t * step) - (lag * step)), ]
      estNames <- row.names(data_lag)
    }

    # Now that we have all future lags, shorten the time-series.
    if (returnMC_full == FALSE & n > t) {
      data_lag <- data_lag[1:((t + 1) * step), ]
      estNames <- row.names(data_lag)[1:((t + 1) * step)]
    }

    # get actual index of Anode
    Anode <- Anode * step + 2 - (lag + 1) * step

    # Get actual index for start, depended on W.
    # (right now, start is set to be the batch, not actual index)
    start <- start * step - (lag + 1) * step + 1
  }

  # TO DO: Later on will need to include more Ws
  # TO DO: This should be parallelized
  outcome <- matrix(nrow = MC, ncol = 1)
  data_lag_origin <- data_lag

  # Need this for the clever covariate.
  res <- matrix(nrow = MC, ncol = 1)

  for (B in 1:MC) {
    data_lag <- data_lag_origin

    # Keep track where we are in the loop
    iter <- 1

    # i always dependent on W iterations. Make everything offset of W
    for (i in seq(start, nrow(data_lag), step)) {

      # Possibly need to skip the update of W and A at the first iteration,
      # depending on from which node we start. Hence the need for iter and ommitting some nodes.

      # Possibly need to condition on earlier time points for calculating our clever covariate.
      # This is our i in the clever covariate.

      if (iter == 1) {
        if (node == "A") {

          # Start at A.
          if (is.null(intervention)) {
            newA <- stats::rbinom(1, 1,
                                  stats::plogis(stats::predict(fit$A,
                                                               data_lag[i + 1, ],
                                                               type = "response")))
            # Update for Y
            data_lag[(i + 2), 1] <- newA
          } else if (!is.null(intervention)) {
            if ((i + 1) == Anode) {
              newA <- stats::rbinom(1, 1, intervention)
              # Update for Y
              data_lag[(i + 2), 1] <- newA
            } else {
              newA <- stats::rbinom(1, 1,
                                    stats::plogis(stats::predict(fit$A,
                                                                 data_lag[i + 1, ],
                                                                 type = "response")))
              # Update for Y
              data_lag[(i + 2), 1] <- newA
            }
          }

          newY <- stats::rbinom(1, 1, stats::plogis(stats::predict(fit$Y,
                                                                   data_lag[i + 2, ],
                                                                   type = "response")))
          # Update for next W
          data_lag[(i + 3), 2] <- newA
          data_lag[(i + 3), 1] <- newY

          iter <- iter + 1
        } else if (node == "Y") {

          # Now we skip both W and A, start generating MC draws from Y.
          newY <- stats::rbinom(1, 1, stats::plogis(stats::predict(fit$Y, data_lag[i + 2, ], type = "response")))
          # Update for next W
          data_lag[(i + 3), 1] <- newY

          iter <- iter + 1
        } else if (node == "W") {

          # Don't skip anything, start with W.
          newW <- stats::rbinom(1, 1, stats::plogis(stats::predict(fit$W, data_lag[i, ], type = "response")))
          # Update for A
          data_lag[(i + 1), 1] <- newW

          if (is.null(intervention)) {
            newA <- stats::rbinom(1, 1, stats::plogis(stats::predict(fit$A, data_lag[i + 1, ], type = "response")))
            # Update for Y
            data_lag[(i + 2), 2] <- newW
            data_lag[(i + 2), 1] <- newA
          } else if (!is.null(intervention)) {
            if ((i + 1) == Anode) {
              newA <- rbinom(1, 1, intervention)
              # Update for Y
              data_lag[(i + 2), 2] <- newW
              data_lag[(i + 2), 1] <- newA
            } else {
              newA <- stats::rbinom(1, 1, stats::plogis(stats::predict(fit$A, data_lag[i + 1, ], type = "response")))
              # Update for Y
              data_lag[(i + 2), 2] <- newW
              data_lag[(i + 2), 1] <- newA
            }
          }

          newY <- stats::rbinom(1, 1, stats::plogis(stats::predict(fit$Y, data_lag[i + 2, ], type = "response")))
          # Update for next W
          # In the last iteration this will be NA... This is ok for now, just one more row.
          data_lag[(i + 3), 2] <- newA
          data_lag[(i + 3), 1] <- newY

          iter <- iter + 1
        }
      } else {

        # For later iterations, do as always.
        # MC draws start from "start" and W node.
        newW <- stats::rbinom(1, 1, stats::plogis(stats::predict(fit$W, data_lag[i, ], type = "response")))
        # Update for A
        data_lag[(i + 1), 1] <- newW
        data_lag[(i + 1), 2] <- newY

        if (is.null(intervention)) {
          newA <- stats::rbinom(1, 1, stats::plogis(stats::predict(fit$A, data_lag[i + 1, ], type = "response")))
          # Update for Y
          data_lag[(i + 2), 2] <- newW
          data_lag[(i + 2), 1] <- newA
        } else if (!is.null(intervention)) {
          if ((i + 1) == Anode) {
            newA <- stats::rbinom(1, 1, intervention)
            # Update for Y
            data_lag[(i + 2), 2] <- newW
            data_lag[(i + 2), 1] <- newA
          } else {
            newA <- stats::rbinom(1, 1, stats::plogis(stats::predict(fit$A, data_lag[i + 1, ], type = "response")))
            # Update for Y
            data_lag[(i + 2), 2] <- newW
            data_lag[(i + 2), 1] <- newA
          }
        }

        newY <- stats::rbinom(1, 1, stats::plogis(stats::predict(fit$Y, data_lag[i + 2, ], type = "response")))
        # Update for next W
        # In the last iteration this will be NA... This is ok for now, just one more row.
        data_lag[(i + 3), 2] <- newA
        data_lag[(i + 3), 1] <- newY

        iter <- iter + 1
      }
    }

    outcome[B, ] <- newY

    if (returnMC == TRUE || returnMC_full == TRUE) {
      data_lag <- data_lag[-1, 1]
      # Collect changed part of the series (from start until the end)
      if (node == "W") {
        retMC[((start + step):nrow(retMC)), B] <- data_lag[(start:length(data_lag))]
        retMC[(1:(start + step - 1)), B] <- data[(1:(start + step - 1)), 1]
      } else if (node == "A") {
        retMC[((start + step + 1):nrow(retMC)), B] <- data_lag[((start + 1):length(data_lag))]
        retMC[(1:(start + step)), B] <- data[(1:(start + step)), 1]
      } else if (node == "Y") {
        retMC[((start + step + 2):nrow(retMC)), B] <- data_lag[((start + 2):length(data_lag))]
        retMC[(1:(start + step + 1)), B] <- data[(1:(start + step + 1)), 1]
      }
    }

    # Get back outcome when specified node at s (or Y^*) was either 1 or 0:
    res[B, ] <- newY
  }

  if (returnMC == TRUE || returnMC_full == TRUE) {
    return(list(
      estimate = mean(outcome), outcome = outcome, intervention = intervention,
      MC = MC, t = t, Anode = Anode,
      s = mean(res, na.rm = TRUE), s_full = res, MCdata = retMC
    ))
  } else {
    return(list(
      estimate = mean(outcome), outcome = outcome, intervention = intervention,
      MC = MC, t = t, Anode = Anode,
      s = mean(res, na.rm = TRUE), s_full = res
    ))
  }
}

