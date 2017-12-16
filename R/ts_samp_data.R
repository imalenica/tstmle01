#' Mock Binary Time Series Data Set
#'
#' A dataset with a simple data structure O = (W, A, Y), where all covariates in
#' the data set are binary (0/1). This is a simple dataset designed specifically
#' to illustrate the TMLE estimation procedure.
#'
#' @format A \code{data.frame} with 3 columns.
#' \describe{
#'   \item{Y}{A binary variable representing an outcome of interest.}
#'   \item{A}{A binary variable representing an intervention of interest.}
#'   \item{W}{A binary variable representing a baseline covariate of interest.}
#' }

"ts_samp_data"

