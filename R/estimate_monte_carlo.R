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
#' @param t Outcome time point of interest. It must be greater than the intervention node A.
#' @param Anode Intervention node.
#' @param intervention Specify %g^*, of %P(A \mid \text{past}).
#' @param lag This is an user impossed dependency lag necessary for the
#'  calculation of the clever covariate in the targeting step. It refers to
#'  %O_i, and it is of the same dimension as the actual C(i), but further in the
#'  past. Default is 0, meaning that it starts right after the node in question.
#'  Also, note that lag should be equal or smaller than start.
#' @param MC How many Monte Carlo samples should be generated.
#' @param returnMC If \code{TRUE}, returns all MC draws up until outcome time.
#' @param returnMC_full If \code{TRUE}, returns all full MC draws.
#' @param clevCov If \code{TRUE}, this MC is used for the calculation of the
#'  clever covariate. Instead of observed data, it used intervened %P^* for further MC draws.
#' @param set Set the s node to either 1 or 0. Used for the clever covariate calculation.
#' @param full \code{TRUE} if full time-series should be used.
#' @param update \code{TRUE}, use updated fits to get Monte Carlo draws.
#'
#' @return An object of class \code{tstmle01}.
#' \describe{
#' \item{estimate}{Mean of the outcome at time t under specified intervention, or no intervention.}
#' \item{outcome}{Outcome at time t for each MC interation.}
#' \item{intervention}{Intervention specified.}
#' \item{MC}{How many Monte Carlo samples should be generated.}
#' \item{Anode}{Intervention node as a function of O_i.}
#' \item{MCdata}{If \code{returnMC} is \code{TRUE}, returns a \code{data.frame} with MC time-series.}
#' \item{set}{Returns what the s node was set to.}
#' }
#'
#' @importFrom stats complete.cases rbinom plogis predict
#' @importFrom Hmisc Lag
#' @importFrom parallel makeCluster detectCores stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach
#'
#' @export
#

mcEst <- function(fit, start = 1, node = "W", t, Anode, intervention = NULL, lag = 0, 
                  MC=100, returnMC = FALSE, returnMC_full = FALSE, clevCov = FALSE, set = NULL, 
                  full=FALSE, update=FALSE) {
  
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
  step <- fit$step
  
  # If estMC is used for clever covariate calculation, use P^*.
  if (clevCov == TRUE && full==FALSE) {
    
    #Already generated O^* with intervention at Anode:
    data <- fit$p_star
    
    #Set s node to "set" (1/0), for Y, A, W.
    #Taking advantage of all MC samples, instead of just averaging over ones that have 1 or 0. 
    #Note: node is the first node that is MCed. So we go back to set the node before.
    if (node == "W") {
      #Set Y:
      #(note, start=s+1 here)
      data[(start * step), ] <- set
    } else if (node == "Y") {
      #Set A:
      #(note, start=s here)
      data[(start * step) + 2, ] <- set
    } else if (node == "A") {
      #Set W:
      #(note, start=s here)
      data[(start * step) + 1, ] <- set
    }
  }else if(clevCov == FALSE && full==TRUE){
    #When simulating p_star and ps for h-densities, use the full data.
    #Shortening should be only for clever cov calculation
    data <- fit$data_orig
  }else if(clevCov == FALSE && full==FALSE){
    data <- fit$data
  }else if(clevCov == TRUE && full==TRUE){
    warning("n for clever covariate and full n might not be the same! Make sure this is the right option.")
  }
  
  # Prepare to return all MCs (up to time t (outcome) generated draws)
  if (returnMC == TRUE & returnMC_full == TRUE) {
    warning("Can't return both a full time-series and a shorter time-series. Returning full time-series.")
  }
  
  # Impose artifical order (C(i)) for the clever covariate calculation.
  # Should default to whatever was estimated in initEst()
  
  if (lag < 0) {
    
    # If lag is negative, move the time series up w.r.t. reference
    
    # This gives further up lag first
    res <- lapply((step * lag + 1):(fit$freqW + step * lag), function(x) {
      Hmisc::Lag(data[, 1], x)
    })
    res <- data.frame(matrix(unlist(res), nrow = length(res[[1]])),
                      stringsAsFactors = FALSE)
    
    data_est_full <- cbind.data.frame(data = data, res)
    data_lag <- data.frame(data_est_full[, -1])
    
    n <- fit$n_true
    
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
    
    #Check the time-series is long enough for this estimation step:
    if(is.na(data_lag[nrow(data_lag),2])){
      stop("The time-series considered is too short for this estimation problem. 
           Consider having a longer time-series, or an earlier outcome point.")
    }
    
    # get actual index of Anode
    Anode <- Anode * step + 2
    
    # Get actual index for start, depended on W.
    # (right now, start is set to be the batch, not actual index)
    start <- start * step + 1
    
    }else if (lag >= 0) {
      
      # This is ok under the assumption orders are the same dimension.
      # Note that "lags" in the clever covariate are defined by the size of the "step".
      # TODO: Change this to an option where we can have varying dimensions for W,A,Y.
      
      res <- lapply((1 + step * lag):(fit$freqW + step * lag), function(x) {
        Hmisc::Lag(data[, 1], x)
      })
      res <- data.frame(matrix(unlist(res), nrow = length(res[[1]])),
                        stringsAsFactors = FALSE)
      
      data_est_full <- cbind.data.frame(data = data, res)
      
      n <- fit$n_true
      
      # Drop time 0 for estimation
      # Notice: drop step*(lag+1) each time
      cc <- stats::complete.cases(data_est_full)
      data_est <- data_est_full[cc, ]
      data_lag <- data.frame(data_est[-1, ])
      data_lag <- data.frame(data_lag[, -1])
      
      # Option for returnMC
      if (clevCov == FALSE & returnMC_full == FALSE) {
        data_lag <- data_lag[1:((t * step) - (lag * step)), ]
        estNames <- row.names(data_lag)
      }
      
      # Now that we have all future lags, shorten the time-series. 
      # Unless we are returning the full MC
      if (returnMC_full == FALSE & n > t) {
        data_lag <- data_lag[1:((t * step)-(lag*step)), ]
        estNames <- row.names(data_lag)[1:((t * step)-(lag*step))]
      }
      
      # get actual index of Anode
      Anode <- Anode * step + 2 - (lag + 1) * step
      
      # Get actual index for start, depended on W.
      # (right now, start is set to be the batch, not actual index)
      start <- start * step - (lag + 1) * step + 1
    }
  
  # TO DO: Later on will need to include more Ws
  
  # i always dependent on W iterations. Make everything offset of W
  #NOTE: last newY is at for data_lag shortened to tau just Y(\tau).
  mainMC<-function(data_lag){
    
    iter <- 1
    
    for (i in seq(start, nrow(data_lag), step)) {
      
      # Possibly need to skip the update of W and A at the first iteration,
      # depending on from which node we start. Hence the need for iter and ommitting some nodes.
      
      # Possibly need to condition on earlier time points for calculating our clever covariate.
      # This is our i in the clever covariate.
      
      if (iter == 1) {
        
        # Start at A.
        if (node == "A") {
          if (is.null(intervention)) {
            # Update for Y
            if(update==TRUE){
              newA <- stats::rbinom(1, 1, stats::plogis(stats::predict(fit$A,data_lag[i + 1, ],type = "response")))
            }else{
              newA <- stats::rbinom(1, 1, stats::predict(fit$A,data_lag[i + 1, ],type = "response")) 
            }
            
            if((i + 2 + (lag*step))>0){
              data_lag[(i + 2 + (lag*step)), 1] <- newA
            }else{
              #Save the obtained newA in place of new.
              data_lag[i + 1, 1]<-newA
            }
          } else if (!is.null(intervention)) {
            if ((i + 1) == Anode) {
              # Update for Y
              newA <- stats::rbinom(1, 1, intervention)
              
              if((i + 2 + (lag*step))>0){
                data_lag[(i + 2 + (lag*step)), 1] <- newA
              }else{
                data_lag[i + 1, 1]<-newA
              }
            } else {
              # Update for Y
              if(update==TRUE){
                newA <- stats::rbinom(1, 1, stats::plogis(stats::predict(fit$A,data_lag[i + 1, ],type = "response")))
              }else{
                newA <- stats::rbinom(1, 1, stats::predict(fit$A,data_lag[i + 1, ],type = "response"))
              }
              
              if((i + 2 + (lag*step))>0){
                data_lag[(i + 2 + (lag*step)), 1] <- newA 
              }else{
                data_lag[i + 1, 1]<-newA
              }
            }
          }
          
          # Update for next W
          if(update==TRUE){
            newY <- stats::rbinom(1, 1, stats::plogis(stats::predict(fit$Y,data_lag[i + 2, ],type = "response")))
          }else{
            newY <- stats::rbinom(1, 1, stats::predict(fit$Y,data_lag[i + 2, ],type = "response"))
          }
          
          if((i + 3 + (lag*step))>0){
            data_lag[(i + 3 + (lag*step)), 2] <- newA
            data_lag[(i + 3 + (lag*step)), 1] <- newY
          }else{
            data_lag[i + 2, 1]<-newY
          }
          
          #Start at Y. 
        }else if (node == "Y") {
          
          # Now we skip both W and A, start generating MC draws from Y.
          # Update for next W
          if(update==TRUE){
            newY <- stats::rbinom(1, 1, stats::plogis(stats::predict(fit$Y, data_lag[i + 2, ], type = "response")))
          }else{
            newY <- stats::rbinom(1, 1, stats::predict(fit$Y, data_lag[i + 2, ], type = "response"))
          }

          if((i + 3 + (lag*step))>0){
            data_lag[(i + 3 + (lag*step)), 1] <- newY
          }else{
            data_lag[i + 2, 1]<-newY
          }
          
          #Start at W.   
        }else if (node == "W") {
          
          # Don't skip anything, start with W.
          # Update for A
          if(update==TRUE){
            newW <- stats::rbinom(1, 1, stats::plogis(stats::predict(fit$W, data_lag[i, ], type = "response")))
          }else{
            newW <- stats::rbinom(1, 1, stats::predict(fit$W, data_lag[i, ], type = "response"))
          }
          
          if((i + 1 + (lag*step))>0){
            data_lag[(i + 1 + (lag*step)), 1] <- newW 
          }else{
            data_lag[i, 1]<-newW
          }
          
          if (is.null(intervention)) {
            # Update for Y
            if(update==TRUE){
              newA <- stats::rbinom(1, 1, stats::plogis(stats::predict(fit$A, data_lag[i + 1, ], type = "response")))
            }else{
              newA <- stats::rbinom(1, 1, stats::predict(fit$A, data_lag[i + 1, ], type = "response"))
            }
            
            if((i + 2 + (lag*step))>0){
              data_lag[(i + 2 + (lag*step)), 2] <- newW
              data_lag[(i + 2 + (lag*step)), 1] <- newA 
            }else{
              data_lag[i + 1, 1]<-newA
            }
          } else if (!is.null(intervention)) {
            if ((i + 1) == Anode) {
              # Update for Y
              newA <- rbinom(1, 1, intervention)
              
              if((i + 2 + (lag*step))>0){
                data_lag[(i + 2 + (lag*step)), 2] <- newW
                data_lag[(i + 2 + (lag*step)), 1] <- newA
              }else{
                data_lag[i + 1, 1]<-newA
              }
            } else {
              # Update for Y
              if(update==TRUE){
                newA <- stats::rbinom(1, 1, stats::plogis(stats::predict(fit$A, data_lag[i + 1, ], type = "response")))
              }else{
                newA <- stats::rbinom(1, 1, stats::predict(fit$A, data_lag[i + 1, ], type = "response"))
              }
              
              if((i + 2 + (lag*step))>0){
                data_lag[(i + 2 + (lag*step)), 2] <- newW
                data_lag[(i + 2 + (lag*step)), 1] <- newA
              }else{
                data_lag[i + 1, 1]<-newA
              }
            }
          }
          
          # Update for next W
          # In the last iteration this will be NA... This is ok for now, just one more row.
          if(update==TRUE){
            newY <- stats::rbinom(1, 1, stats::plogis(stats::predict(fit$Y, data_lag[i + 2, ], type = "response")))
          }else{
            newY <- stats::rbinom(1, 1, stats::predict(fit$Y, data_lag[i + 2, ], type = "response"))
          }

          if((i + 3 + (lag*step))>0){
            data_lag[(i + 3 + (lag*step)), 2] <- newA
            data_lag[(i + 3 + (lag*step)), 1] <- newY 
          }else{
            data_lag[i + 2, 1]<-newY
          }
          
        }
        
        #Later on iterations.  
      }else {
        
        # For later iterations, do as always.
        # MC draws start from "start" and W node.
        
        # Update for A
        if(update==TRUE){
          newW <- stats::rbinom(1, 1, stats:plogis(stats::predict(fit$W, data_lag[i, ], type = "response")))
        }else{
          newW <- stats::rbinom(1, 1, stats::predict(fit$W, data_lag[i, ], type = "response"))
        }

        if((i + 1 + (lag*step))>0){
          data_lag[(i + 1 + (lag*step)), 1] <- newW
          data_lag[(i + 1 + (lag*step)), 2] <- newY 
        }else{
          data_lag[i, 1]<-newW
        }
        
        #Intervention=NULL
        if (is.null(intervention)) {
          # Update for Y
          if(update==TRUE){
            newA <- stats::rbinom(1, 1, stats::plogis(stats::predict(fit$A, data_lag[i + 1, ], type = "response")))
          }else{
            newA <- stats::rbinom(1, 1, stats::predict(fit$A, data_lag[i + 1, ], type = "response"))
          }

          if((i + 2 + (lag*step))>0){
            data_lag[(i + 2 + (lag*step)), 2] <- newW
            data_lag[(i + 2 + (lag*step)), 1] <- newA
          }else{
            data_lag[i + 1, 1]<-newA
          }
        } else if (!is.null(intervention)) {
          
          #Intervention at the current node batch
          if ((i + 1) == Anode) {
            # Update for Y
            newA <- stats::rbinom(1, 1, intervention)
            
            if((i + 2 + (lag*step))>0){
              data_lag[(i + 2 + (lag*step)), 2] <- newW
              data_lag[(i + 2 + (lag*step)), 1] <- newA
            }else{
              data_lag[i + 1, 1]<-newA
            }
            
            #There is an intervention but not in this batch  
          } else {
            # Update for Y
            if(update==TRUE){
              newA <- stats::rbinom(1, 1, stats::plogis(stats::predict(fit$A, data_lag[i + 1, ], type = "response")))
            }else{
              newA <- stats::rbinom(1, 1, stats::predict(fit$A, data_lag[i + 1, ], type = "response"))
            }

            if((i + 2 + (lag*step))>0){
              data_lag[(i + 2 + (lag*step)), 2] <- newW
              data_lag[(i + 2 + (lag*step)), 1] <- newA
            }else{
              data_lag[i + 1, 1]<-newA
            }
          }
        }
        
        # Update for next W
        # In the last iteration this will be NA... This is ok for now, just one more row.
        if(update==TRUE){
          newY <- stats::rbinom(1, 1, stats::plogis(stats::predict(fit$Y, data_lag[i + 2, ], type = "response")))
        }else{
          newY <- stats::rbinom(1, 1, stats::predict(fit$Y, data_lag[i + 2, ], type = "response"))
        }

        if((i + 3 + (lag*step))>0){
          data_lag[(i + 3 + (lag*step)), 2] <- newA
          data_lag[(i + 3 + (lag*step)), 1] <- newY 
        }else{
          data_lag[i + 2, 1]<-newY
        }
        
      }
      
      iter <- iter + 1
    }
    
    if(returnMC == TRUE || returnMC_full == TRUE){
      
      if(lag>=0){
        if(node=="A"){
          data_lag_MC <- data.frame(data=data_lag[(lag*step+3):nrow(data_lag), 1])
          row.names(data_lag_MC)<-row.names(data_lag)[2:(nrow(data_lag)-(lag*step+1))]
          
          #data_lag is truncated to go from start. Add what's missing
          missData<-data.frame(data=data[1:(step+1+(abs(lag)*step)),])
          row.names(missData)<-row.names(data)[1:(step+1+(abs(lag)*step))]
          data_lag_MC<-rbind.data.frame(missData,data_lag_MC)
        }else if(node=="W"){
          data_lag_MC <- data.frame(data=data_lag[(lag*step+2):nrow(data_lag), 1])
          row.names(data_lag_MC)<-row.names(data_lag)[1:(nrow(data_lag)-(lag*step+1))]
          
          #data_lag is truncated to go from start. Add what's missing
          missData<-data.frame(data=data[1:(step+(abs(lag)*step)),])
          row.names(missData)<-row.names(data)[1:(step+(abs(lag)*step))]
          data_lag_MC<-rbind.data.frame(missData,data_lag_MC)
        }else if(node=="Y"){
          data_lag_MC <- data.frame(data=data_lag[(lag*step+4):nrow(data_lag), 1])
          row.names(data_lag_MC)<-row.names(data_lag)[3:(nrow(data_lag)-(lag*step+1))]
          
          #data_lag is truncated to go from start. Add what's missing
          missData<-data.frame(data=data[1:(step+2+(abs(lag)*step)),])
          row.names(missData)<-row.names(data)[1:(step+2+(abs(lag)*step))]
          data_lag_MC<-rbind.data.frame(missData,data_lag_MC)
        }
        #TO DO: Investigate for starting node W,A,Y  
      }else if(lag<0){
        
        if((nrow(data_lag)+(lag*step+1))>0){
          #Goes until the end
          data_lag_MC<-data.frame(data=data_lag[1:(nrow(data_lag)+(lag*step+1)),1])
          row.names(data_lag_MC)<-row.names(data_lag)[-(lag*step):nrow(data_lag)]
          
          #data_lag is truncated to go from start. Add what's missing
          missData<-data.frame(data=data[1:-(lag*step+1),])
          row.names(missData)<-row.names(data)[1:-(lag*step+1)]
          data_lag_MC<-rbind.data.frame(missData,data_lag_MC)
        }else{
          #No updating occured: return what was calculated based on the original ts
          data_lag_MC<-data.frame(data=data_lag[start:nrow(data_lag),1])
          data_lag_MC<-rbind.data.frame(data.frame(data=data_est_full[1:(start-1),1]),data_lag_MC)
          row.names(data_lag_MC)<-row.names(data_lag)
        }
        
      }
      
      return(list(newY=newY,dataMC=data_lag_MC))
      
    }else {
      return(list(newY=newY))
    }
    
  }
  
  #Set up to run in parallel
  no_cores <- parallel::detectCores() - 1
  myCluster <- parallel::makeCluster(no_cores)
  doParallel::registerDoParallel(myCluster)
  
  #Tracks Y(\tau) and MC data
  outcome_all<-foreach::foreach(1:MC, .combine = c) %dopar% mainMC(data_lag=data_lag)
  base::on.exit(parallel::stopCluster(myCluster))
  
  if(returnMC == TRUE){
    
    a <- 1:(2*MC)
    
    #Save Y(\tau)
    outcome<-unlist(lapply(a[seq(1, length(a), 2)], function(x) outcome_all[x]$newY))
    
    #Save MC data
    dataMC<-lapply(a[seq(2, length(a), 2)], function(x) outcome_all[x]$dataMC)
    dataMC <- data.frame(matrix(unlist(dataMC), nrow=(step*(t+1)), byrow=F))
    names(dataMC)<-paste0("MC",1:MC)
    row.names(dataMC)<-row.names(data)[1:(step*(t+1))]
    
  }else if(returnMC_full == TRUE){
    
    a <- 1:(2*MC)
    
    #Save Y(\tau)
    outcome<-unlist(lapply(a[seq(1, length(a), 2)], function(x) outcome_all[x]$newY))
    
    #Save MC data
    dataMC<-lapply(a[seq(2, length(a), 2)], function(x) outcome_all[x]$dataMC)
    dataMC <- data.frame(matrix(unlist(dataMC), nrow=(step*(n+1)), byrow=F))
    names(dataMC)<-paste0("MC",1:MC)
    row.names(dataMC)<-row.names(data)
    
  }else{
    #Save Y(\tau)
    outcome<-unlist(lapply(1:MC, function(x) outcome_all[x]$newY))
  }
  
  if (returnMC == TRUE || returnMC_full == TRUE) {
    return(list(estimate = mean(outcome), outcome = outcome, intervention = intervention,
                MC = MC, t = t, Anode = Anode, set=set, MCdata = dataMC))
    
  }else {
    return(list(estimate = mean(outcome), outcome = outcome, intervention = intervention,
                MC = MC, t = t, Anode = Anode, set=set))
  }
  
}