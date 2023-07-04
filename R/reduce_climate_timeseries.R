# Functions in this file do basic smoothing, integration, and other 
# manipulations high time-resolution data

#' Fit a spline to a collection of timeseries, with a possible roughness penalty
#'
#' @param data_to_smooth A data frame. Rows correspond to different time points
#' and columns correspond to different variables.
#' @param timepoints A vector (numeric). Corresponds to the measurement 
#' timepoints of each variable.
#' @param nbasis_in Integer. The number of basis functions to use.
#' @param lambda_in Float. Defaults to 0. The value of the smoothing parameter,
#' lambda.
#' @param basis_type String. Defaults to "fourier". Can be either "fourier" or
#' "bspline" to select different families of basis functions.
#' @param norder_in Integer. Defaults to 4. Only relevant if basis_type = 
#' TRUE. Controls the order of the bspline functions.
#'
#' @return An fdSmooth object.
fits_from_params <- function(data_to_smooth,
                             timepoints,
                             nbasis_in,
                             lambda_in = 0,
                             basis_type = "fourier",
                             norder_in = 4) {
  range_of_timeseries <- c(min(timepoints), max(timepoints))
  
  penalty_basis <- switch(
    basis_type,
    fourier = fda::create.fourier.basis(
      rangeval = range_of_timeseries,
      nbasis = nbasis_in,
      
    ),
    bspline = fda::create.bspline.basis(
      rangeval = range_of_timeseries,
      nbasis = nbasis_in,
      norder = norder_in
    )
  )
  
  if (is.null(penalty_basis)) {
    stop('basis_type not recognised - must be "fourier" or "bspline"')
  }
  
  
  # Define a linear differential operator
  second_deriv_operator = fda::int2Lfd(2)
  
  fd_params = fda::fdPar(penalty_basis, second_deriv_operator, lambda_in)
  fd_outputs = fda::smooth.basis(
    timepoints,
    data_to_smooth, 
    fd_params
  )
  
  return(fd_outputs)
  
}
