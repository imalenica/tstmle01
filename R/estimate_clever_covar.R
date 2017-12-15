#' Clever Covariate calculation
#'
#' This function calculates the clever covariate for each time point.
#'
#' @param fit \code{fit} object obtained by \code{initEst}.
#' @param t Outcome time point of interest. It must be greater than the intervention node A.
#' @param Anode Intervention node.
#' @param intervention Specify %g^*, of %P(A \mid \text{past}). Right now, this supports only 1/0 type interventions.
#' @param MC How many Monte Carlo samples should be generated.
#' @param B How many samples to draw from P, as part of the h-density estimation.
#' @param N How many sample to draw from %P^*, as part of the h-density estimation.
#'
#' @return An object of class \code{tstmle}.
#' \describe{
#' \item{Hy}{Clever covariate for Y component of the likelihood.}
#' \item{Ha}{Clever covariate for A component of the likelihood.}
#' \item{Hw}{Clever covariate for W component of the likelihood.}
#' \item{Dbar}{Efficient Influence Curve.}
#' \item{h_star}{Empirical estimate of h-density under the intervention g^*.}
#' \item{h}{Empirical estimate of h-density under no intervention.}
#' }
#'
#' @export
#

cleverCov <- function(fit, t, Anode, intervention = 1, B = 100, N = 100, MC = 100) {
  
  step <- length(grep("_1$", row.names(fit$data), value = TRUE))
  n <- nrow(fit$data) / step - 1

  # Generate clever covariates for each of the likelihood components: W, A, Y
  # and efficient influence function (EIF)
  Hy_cc <- matrix(nrow = n, ncol = 1)
  Ha_cc <- matrix(nrow = n, ncol = 1)
  Hw_cc <- matrix(nrow = n, ncol = 1)

  # h-density ratio components
  hy_res <- matrix(nrow = n, ncol = 1)
  ha_res <- matrix(nrow = n, ncol = 1)
  hw_res <- matrix(nrow = n, ncol = 1)

  hy_star_res <- matrix(nrow = n, ncol = 1)
  ha_star_res <- matrix(nrow = n, ncol = 1)
  hw_star_res <- matrix(nrow = n, ncol = 1)

  # EIC:
  D <- matrix(nrow = n, ncol = 1)

  # TO DO: Probably need some kind of an internal seed for these computations.
  # Generate our P^*, intervening only on Anode.
  p_star <- mcEst(fit, start = Anode, node = "A", t = t, Anode = Anode,
                  intervention = intervention, MC = 1, returnMC_full = TRUE)

  # Add intervened data to fit:
  fit[["p_star"]] <- p_star$MCdata

  for (i in seq_len(n)) {
    
    Hy_diff <- 0
    Ha_diff <- 0
    Hw_diff <- 0

    hy_star <- 0
    ha_star <- 0
    hw_star <- 0

    for (s in seq_len(t)) {

      if (s < i) {
        # Look into the future to condition on.
        if (s == t) {
          #Y; If s=t, then just return Y^* for Hy.
          Hy_diff_add <- p_star$MCdata[length(p_star$MCdata), ]
          
          #A
          Ha1 <- mcEst(fit, start = s, node = "Y", t = t, Anode = Anode, 
                       lag = (-1 * (i - s)), MC = MC, clevCov = TRUE, set = 1)
          Ha0 <- mcEst(fit, start = s, node = "Y", t = t, Anode = Anode,
                       lag = (-1 * (i - s)), MC = MC, clevCov = TRUE, set = 0)
          Ha_diff_add <- Ha1$s - Ha0$s

          #W
          Hw1 <- mcEst(fit, start = s, node = "A", t = t, Anode = Anode,
                       lag = (-1 * (i - s)), MC = MC, clevCov = TRUE, set = 1)
          Hw0 <- mcEst(fit, start = s, node = "A", t = t, Anode = Anode,
                       lag = (-1 * (i - s)), MC = MC, clevCov = TRUE, set = 0)
          Hw_diff_add <- Hw1$s - Hw0$s
          
        } else {
          #Y
          Hy1 <- mcEst(fit, start = s + 1, node = "W", t = t, Anode = Anode,
                       lag = (-1 * (i - s)), MC = MC, clevCov = TRUE, set = 1)
          Hy0 <- mcEst(fit, start = s + 1, node = "W", t = t, Anode = Anode,
                       lag = (-1 * (i - s)), MC = MC, clevCov = TRUE, set = 0)
          Hy_diff_add <- Hy1$s - Hy0$s
          
          #A
          Ha1 <- mcEst(fit, start = s, node = "Y", t = t, Anode = Anode,
                       lag = (-1 * (i - s)), MC = MC, clevCov = TRUE, set = 1)
          Ha0 <- mcEst(fit, start = s, node = "Y", t = t, Anode = Anode,
                       lag = (-1 * (i - s)), MC = MC, clevCov = TRUE, set = 0)
          Ha_diff_add <- Ha1$s - Ha0$s

          #W
          Hw1 <- mcEst(fit, start = s, node = "A", t = t, Anode = Anode,
                       lag = (-1 * (i - s)), MC = MC, clevCov = TRUE, set = 1)
          Hw0 <- mcEst(fit, start = s, node = "A", t = t, Anode = Anode,
                       lag = (-1 * (i - s)), MC = MC, clevCov = TRUE, set = 0)
          Hw_diff_add <- Hw1$s - Hw0$s
        }
        
      } else if (s == i) {
        # Usual setting.
        if (s == t) {
          
          #Y; If s=t, then just return Y^* for Hy.
          Hy_diff_add <- p_star$MCdata[length(p_star$MCdata), ]

          #A
          Ha1 <- mcEst(fit, start = s, node = "Y", t = t, Anode = Anode,
                       lag = 0, MC = MC, clevCov = TRUE, set = 1)
          Ha0 <- mcEst(fit, start = s, node = "Y", t = t, Anode = Anode,
                       lag = 0, MC = MC, clevCov = TRUE, set = 0)
          Ha_diff_add <- Ha1$s - Ha0$s

          #W
          Hw1 <- mcEst(fit, start = s, node = "A", t = t, Anode = Anode,
                       lag = 0, MC = MC, clevCov = TRUE, set = 1)
          Hw0 <- mcEst(fit, start = s, node = "A", t = t, Anode = Anode,
                       lag = 0, MC = MC, clevCov = TRUE, set = 0)
          Hw_diff_add <- Hw1$s - Hw0$s
          
        } else {
          #Y
          Hy1 <- mcEst(fit, start = s + 1, node = "W", t = t, Anode = Anode,
                       lag = 0, MC = MC, clevCov = TRUE, set = 1)
          Hy0 <- mcEst(fit, start = s + 1, node = "W", t = t, Anode = Anode,
                       lag = 0, MC = MC, clevCov = TRUE, set = 0)
          Hy_diff_add <- Hy1$s - Hy0$s

          #A
          Ha1 <- mcEst(fit, start = s, node = "Y", t = t, Anode = Anode,
                       lag = 0, MC = MC, clevCov = TRUE, set = 1)
          Ha0 <- mcEst(fit, start = s, node = "Y", t = t, Anode = Anode,
                       lag = 0, MC = MC, clevCov = TRUE, set = 0)
          Ha_diff_add <- Ha1$s - Ha0$s

          #W
          Hw1 <- mcEst(fit, start = s, node = "A", t = t, Anode = Anode,
                       lag = 0, MC = MC, clevCov = TRUE, set = 1)
          Hw0 <- mcEst(fit, start = s, node = "A", t = t, Anode = Anode,
                       lag = 0, MC = MC, clevCov = TRUE, set = 0)
          Hw_diff_add <- Hw1$s - Hw0$s
        }
        
      } else if (s > i) {
        # Look further in the past since s>i.
        if (s == t) {
          
          #Y; If s=t, then just return Y^* for Hy.
          Hy_diff_add <- p_star$MCdata[length(p_star$MCdata), ]

          #A
          Ha1 <- mcEst(fit, start = s, node = "Y", t = t, Anode = Anode,
                       lag = (s - i), MC = MC, clevCov = TRUE, set = 1)
          Ha0 <- mcEst(fit, start = s, node = "Y", t = t, Anode = Anode,
                       lag = (s - i), MC = MC, clevCov = TRUE, set = 0)
          Ha_diff_add <- Ha1$s - Ha0$s

          #W
          Hw1 <- mcEst(fit, start = s, node = "A", t = t, Anode = Anode,
                       lag = (s - i), MC = MC, clevCov = TRUE, set = 1)
          Hw0 <- mcEst(fit, start = s, node = "A", t = t, Anode = Anode,
                       lag = (s - i), MC = MC, clevCov = TRUE, set = 0)
          Hw_diff_add <- Hw1$s - Hw0$s
          
        } else {
          #Y
          Hy1 <- mcEst(fit, start = s + 1, node = "W", t = t, Anode = Anode,
                       lag = (s - i), MC = MC, clevCov = TRUE, set = 1)
          Hy0 <- mcEst(fit, start = s + 1, node = "W", t = t, Anode = Anode,
                       lag = (s - i), MC = MC, clevCov = TRUE, set = 0)
          Hy_diff_add <- Hy1$s - Hy0$s

          #A
          Ha1 <- mcEst(fit, start = s, node = "Y", t = t, Anode = Anode,
                       lag = (s - i), MC = MC, clevCov = TRUE, set = 1)
          Ha0 <- mcEst(fit, start = s, node = "Y", t = t, Anode = Anode,
                       lag = (s - i), MC = MC, clevCov = TRUE, set = 0)
          Ha_diff_add <- Ha1$s - Ha0$s

          #W
          Hw1 <- mcEst(fit, start = s, node = "A", t = t, Anode = Anode,
                       lag = (s - i), MC = MC, clevCov = TRUE, set = 1)
          Hw0 <- mcEst(fit, start = s, node = "A", t = t, Anode = Anode,
                       lag = (s - i), MC = MC, clevCov = TRUE, set = 0)
          Hw_diff_add <- Hw1$s - Hw0$s
        }
      }

      Hy_diff <- Hy_diff + Hy_diff_add
      Ha_diff <- Ha_diff + Ha_diff_add
      Hw_diff <- Hw_diff + Hw_diff_add

      # Get h^*
      h_star_est <- hstarEst(fit, s, i, B, t, Anode, intervention)

      # Sum over all s:
      hy_star <- hy_star + h_star_est$h_cy_star
      ha_star <- ha_star + h_star_est$h_ca_star
      hw_star <- hw_star + h_star_est$h_cw_star
    }

    # Get h
    h_init <- hEst(fit, i, B = N, t)

    hy_res[i, ] <- h_init$h_cy
    ha_res[i, ] <- h_init$h_ca
    hw_res[i, ] <- h_init$h_cw

    # Save sums for h^*
    hy_star_res[i, ] <- hy_star
    ha_star_res[i, ] <- ha_star
    hw_star_res[i, ] <- hw_star

    # Notice that this will give a LOT of weight to early time-points if t<N
    # (whenever there is a chance i=s)

    # Generate ratio for each component of the likelihood for "sample" time i:
    hy_ratio <- hy_star_res[i, ] / hy_res[i, ]
    ha_ratio <- ha_star_res[i, ] / ha_res[i, ]
    hw_ratio <- hw_star_res[i, ] / hw_res[i, ]

    # Store clever covariates for each time point
    Hy_cc[i, ] <- Hy_diff * hy_ratio
    
    # On intervention node 0:
    if (i == Anode) {
      Ha_cc[i, ] <- 0
    } else {
      Ha_cc[i, ] <- Ha_diff * ha_ratio
    }
    
    Hw_cc[i, ] <- Hw_diff * hw_ratio

    # Calculate the EIC:
    preds <- getPred(fit, i)
    D[i, ] <- Hy_cc[i, ] * (preds$Y - preds$Y_pred) + Ha_cc[i, ] *
      (preds$A - preds$A_pred) + Hw_cc[i, ] * (preds$W - preds$W_pred)
  }

  return(list(Hy = Hy_cc, Ha = Ha_cc, Hw = Hw_cc, Dbar = D,
              h_star = list(Hy_star = hy_star_res, Ha_star = ha_star_res,
                            Hw_star = hw_star_res),
              h = list(Hy = hy_res, Ha = ha_res, Hw = hw_res)
        ))
}

