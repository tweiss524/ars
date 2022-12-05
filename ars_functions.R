

check_input <- function(f, f_prime, n, D_bounds, k, x_init) {
  
  assertthat::assert_that(is.function(f),
                          msg = "f must be a function")
  assertthat::assert_that(is.function(f_prime),
                          msg = "f_prime must be a function")
  assertthat::assert_that(is.integer(n),
                          msg = "n must be an integer")
  assertthat::assert_that(is.vector(D_bounds) && (length(D_bounds) == 2) &&
              (is.numeric(D_bounds)),
            msg = "Bounds must be numieric vector of length 2")
  
  if (D_bounds[1] > D_bounds[2]) {
    
    D_bounds <- sort(D_bounds)
    
    warning(paste0("Lower bound must be smaller than upper bound. ", 
            sprintf("New bounds are (%s, %s)", D_bounds[1], D_bounds[2])))
    
  }
  
  if (D_bounds[1] == D_bounds[2]) {
    
    D_bounds <- c(-Inf, Inf)
    
    warning(paste0("Lower bound must be smaller than upper bound. ", 
                   sprintf("New bounds are (%s, %s)", D_bounds[1], D_bounds[2])))
    
  }
  
  assertthat::assert_that((is.integer(k)) && (k > 0), msg = "k must be a positive integer")
  
  assertthat::assert_that((x_init[1] <= D_bounds[1]) && (x_init[2] <= D_bounds[2]),
             msg = "Inital point must be inside bounds")
  
}

################################################################################

h <- function(x) {
  return(log(f(x)))
}

################################################################################


hprime <- function(x) {
  der <- numDeriv::grad(h, x)
  return(der)
}

###############################################################################


check_log_concave <- function(x){
  ep <- 1e-8
  n <- length(x)
  assertthat::assert_that(x[2:n] - x[1:n-1] < ep,
                          msg = "Function is not log-concave")
}


###############################################################################







