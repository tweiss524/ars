

check_input <- function(f, n, bounds, k, x_init, mu) {
  
  assertthat::assert_that(is.function(f), msg = "f must be a function")
  assertthat::assert_that(is.numeric(n), msg = "n must be an integer")
  assertthat::assert_that(is.numeric(x_init), msg = "x_init must be an integer")
  assertthat::assert_that(is.numeric(mu), msg = "mu must be numeric")
  assertthat::assert_that(is.vector(bounds) && (length(bounds) == 2) && (is.numeric(bounds)), msg = "Bounds must be numeric vector of length 2")
  
  if (bounds[1] > bounds[2]) {
    
    bounds <<- sort(bounds)
    bounds <- sort(bounds)
    
    warning(paste0("Lower bound must be smaller than upper bound. ", 
            sprintf("New bounds are (%s, %s)", bounds[1], bounds[2])))
    
  }
  
  if (bounds[1] == bounds[2]) {
    
    bounds <<- c(-Inf, Inf)
    bounds <- c(-Inf, Inf)
    
    warning(paste0("Lower bound must be smaller than upper bound. ", 
                   sprintf("New bounds are (%s, %s)", bounds[1], bounds[2])))
    
  }
  
  assertthat::assert_that((is.numeric(k)) && (k > 0), msg = "k must be a positive integer")
  
  assertthat::assert_that((x_init >= bounds[1]) && (x_init <= bounds[2]), msg = "Inital point must be inside bounds")
  
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
  assertthat::assert_that(x[2:n] - x[1:n-1] < ep, msg = "Function is not log-concave")
}


###############################################################################

z <- function(x1, x2, h, hprime) {
  
  return((h(x2) - h(x1) - x2 * hprime(x2) + x1* hprime(x1))/(hprime(x1)-hprime(x2)))
  
}


initialize_Tk <- function(k, bounds, hprime, mu) {
  
  if((bounds[1] == -Inf) && (bounds[2] == Inf)) {
    x <- mu - 10
    y <- mu + 10
    
    while((!is.finite(hprime(x))) && (!is.finite(hprime(y)))) {
      x <- x + 1
      y <- y - 1
      
      if (x >= y) {
        stop("Invalid Bounds.")
      }
      
    }
    
    bounds[1] <- x
    bounds[2] <- y
    
  }
  
  if((bounds[1] == -Inf) && (bounds[2] != Inf)) {
    
    x <- bounds[2] - 20
    
    while((hprime(x) <= 0) && (hprime(x) != Inf) && (x < bounds[2])) {
      x <- x + 1
    }
    
    if((hprime(x) == Inf) | (x >= bounds[2])) {
      stop("Invalid Bounds.")
    }
    
    bounds[1] <- x
    
  }
  
  if((bounds[1] != -Inf) && (bounds[2] == Inf)) {
    
    y <- bounds[1] + 20
    
    while((hprime(y) >= 0) && (hprime(y) != -Inf) && (y > bounds[1])) {
      y <- y - 1
    }
    
    if((hprime(y) == -Inf) | (y <= bounds[1])) {
      stop("Invalid Bounds.")
    }
    
    bounds[2] <- y
    
  }
  
  if((bounds[1] != -Inf) && (bounds[2] != Inf)) {
    
    if((!is.finite(hprime(bounds[1]))) || (!is.finite(hprime(bounds[2])))) {
      stop("Invalid bounds.")
    }
  }
  
  Tk <- seq(bounds[1], bounds[2], length.out = k)
}


