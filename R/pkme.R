#' Solve PK or PBPK model in matrix form
#'
#' Takes in initial condition matrix (single column vector of initial concentrations) and a rate
#' coefficient matrix (flow/volume), dose profile (number of doses, start and stop time of each dose)
#' along with any parameters needed to solve the system, and outputs a matrix of the concentrations
#' in each compartment (rows) versus the time course (columns) for a given step size "h". Each column
#' represents a predicted data point at every "2*h"
#' @param minute The number of minutes the model will simulate out
#' @param h The step size specified for Simpson's 1/3 rule
#' @param rate_coeff A matrix of rates - flow/volume corresponding to the equation/compartment (row) and concentration to be multiplied by (column)
#' @param init_condition A column vector of initial conditions (concentrations) for each equation/compartment (row)
#' @param pars A matrix of parameters defining the variables in the rate coefficient matrix and the initial conditions, if necessary
#' @param number_of_doses A scalar indicating number of doses
#' @param dose A vector of doses corresponding to the "ith" start and stop time vector
#' @param start_time A vector of start times for an "i" number of doses
#' @param stop_time A vector of stop times for an "i" number of doses
#' @param dose_comp Designate which compartment to place the dose -> default is 1
#' @return x_matrix_save A matrix of the concentrations in each compartment (rows) versus the time course (columns) for a given step size "h". Each column represents a predicted data point at every "2*h"
#' @export


pkme <- function(minute, h, rate_coeff, init_condition, pars, number_of_doses, dose, start_time, stop_time, dose_comp = 1) {
  t_index_max = (minute / h) + 1
  tol_value = 10^-500
  matrix_size = nrow(rate_coeff)
  injection_rate <- matrix(0, nrow = t_index_max, ncol=1)
  for (t_index in 1:t_index_max) {
   t <- (t_index - 1) * h
   for (i in 1:number_of_doses) {
     if (t <= stop_time[i] && t > start_time[i]) {injection_rate[t_index,] <- dose[i] / (stop_time[i] - start_time[i])}
   }
  }


  if (nrow(rate_coeff) != ncol(rate_coeff)) {
    stop('matrix must be square')
  }

  if (nrow(rate_coeff) != nrow(init_condition)) {
    stop('initial conditions must be same size as rate coefficient matrix')
  }
  # eigenvalue/eigenvector solution

  # eigenvalues and eigenvectors of coefficient matrix
  eigen_matrix <- eigen(rate_coeff)

  # print eigenvalues and eigenvectors
  eigen_matrix$values
  eigen_matrix$vectors

  # separate real and imaginary parts of eigenvectors
  a_real_eigenvec <- Re(eigen_matrix$vectors)
  b_imag_eigenvec <- Im(eigen_matrix$vectors)

  # print real and complex eigenvectors
  a_real_eigenvec
  b_imag_eigenvec

  # separate real and imaginary parts of eigenvalues
  lambda <- Re(eigen_matrix$values)
  mu <- Im(eigen_matrix$values)

  # print real and complex eigenvalues
  lambda
  mu

  # initialize matrices
  conc_v_time <- vector()
  u_matrix <- matrix(0, nrow = matrix_size, ncol = 1)
  g_vector <- matrix(0, nrow = matrix_size, ncol = 1)
  u_prime_matrix <- matrix()
  u_prime_matrix_save <- matrix(0, nrow = matrix_size, ncol = 1)
  fund_matrix <- matrix(0, nrow = matrix_size, ncol = matrix_size)

  t = 0
  ## calculate fundamental matrix and inverse of fundamental matrix for 0, h, 2*h
  while (t <= h*2) {
    i <- 1
    while (i <= matrix_size){
      j <- 1
      while (j <= matrix_size) {
        if (mu[j] == 0) {
          fund_matrix[i,j] <- exp(lambda[j]*t)*(a_real_eigenvec[i,j] * cos(mu[j]*t))
          j <- j + 1
        }
        else {
          fund_matrix[i,j] <- exp(lambda[j]*t) * (a_real_eigenvec[i,j] * cos(mu[j]*t) - b_imag_eigenvec[i,j] * sin(mu[j]*t))
          j <- j + 1
          fund_matrix[i,j] <- exp(lambda[j]*t) * (a_real_eigenvec[i,j] * sin(mu[j]*t) + b_imag_eigenvec[i,j] * cos(mu[j]*t))
          j <- j + 1
        }
      }
      i = i + 1
    }

    if (t == 0) {
      t0_fund_matrix <- fund_matrix
    } else if (t == h * 1) {
      t1_fund_matrix <- fund_matrix
    } else if (t == h * 2 ) {
      t2_fund_matrix <- fund_matrix
    } else {
      fund_matrix <- fund_matrix
    }

    t = t + h
  }

  t0_inv_fund_matrix = solve(t0_fund_matrix, tol = tol_value)
  t1_inv_fund_matrix = solve(t1_fund_matrix, tol = tol_value)
  t2_inv_fund_matrix = solve(t2_fund_matrix, tol = tol_value)

  ## function starts here

  for (t_index in seq(from = 1, to = t_index_max, by = 1)) {
    g_vector[dose_comp,1] <- injection_rate[t_index]

    u_prime_matrix0 <- t0_inv_fund_matrix %*% g_vector
    u_prime_matrix1 <- t1_inv_fund_matrix %*% g_vector
    u_prime_matrix2 <- t2_inv_fund_matrix %*% g_vector

    ## initial time,
    ## run once to get constant vector
    t0_inv_fund_matrix = solve(t0_fund_matrix, tol = tol_value)
    # constants for a given initial condition at t = 0
    int_constants =   t0_inv_fund_matrix %*% init_condition

    ## Simpson's 1/3 rule for 0 to 0.2
    sim_soln <- (h/3) * (u_prime_matrix0 + (4 * u_prime_matrix1) + u_prime_matrix2)
    sim_soln <- matrix(sim_soln)

    ## u_matrix should be term by term u_prime_matrix integrated w.r.t. t and bounds 0 to time point
    x_matrix = fund_matrix %*% int_constants + fund_matrix %*% sim_soln
    conc_v_time <- cbind(conc_v_time, x_matrix)

    ## x_matrix solution becomes initial condition for next tau
    init_condition = x_matrix
   }
   return(conc_v_time)
  }
