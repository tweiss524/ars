
ars <- function(f, n = 1000,
                bounds = c(-Inf, Inf), k = 20,
                x_init = 1) {
  
  
  check_input <- function(f, n, bounds, k) {
    
    assertthat::assert_that(is.function(f), msg = "f must be a function")
    assertthat::assert_that(is.numeric(n), msg = "n must be an integer")
    #assertthat::assert_that(is.numeric(x_init), msg = "x_init must be an integer")
  
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
    
    #assertthat::assert_that((x_init >= bounds[1]) && (x_init <= bounds[2]), msg = "Inital point must be inside bounds")
    
  }
  
  ################################################################################
  
  check_log_concave <- function(x){
    ep <- 1e-8
    n <- length(x)
    assertthat::assert_that(sum(x[2:n] - x[1:n-1] < ep) == n-1, msg = "Function is not log-concave")
  }
  
  
  ###############################################################################
  
  calc_z <- function(x, h_tk, hprime_tk) {
    n <- length(x)
    x1 <- x[1:(n-1)]
    x2 <- x[2:n]
    
    return((h_tk[2:n] - h_tk[1:(n-1)] - x[2:n] * hprime_tk[2:n] + x[1:(n-1)]* hprime_tk[1:(n-1)])/(hprime_tk[1:(n-1)]-hprime_tk[2:n]))
    
  }
  
  ########################################################################
  
  initialize_abcissae <- function(k, bounds, hprime) {
    
    if((bounds[1] == -Inf) && (bounds[2] == Inf)) {
      x <- mu - 10
      y <- mu + 10
      
      while((!is.finite(hprime(x))) && (!is.finite(hprime(y)))) {
        x <- x + 1
        y <- y - 1
        
        if (x >= y) {
          stop("Invalid Bounds case (-Inf, Inf).")
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
        stop("Invalid Bounds case (-Inf, x).")
      }
      
      bounds[1] <- x
      
    }
    
    if((bounds[1] != -Inf) && (bounds[2] == Inf)) {
      
      y <- bounds[1] + 20
      
      while((hprime(y) >= 0) && (hprime(y) != -Inf) && (y > bounds[1])) {
        y <- y - 1
      }
      
      if((hprime(y) == -Inf) | (y <= bounds[1])) {
        stop("Invalid Bounds case (x, Inf).")
      }
      
      bounds[2] <- y
      
    }
    
    if((bounds[1] != -Inf) && (bounds[2] != Inf)) {
      
      if((!is.finite(hprime(bounds[1]))) || (!is.finite(hprime(bounds[2])))) {
        stop("Function f not defined on bounds case (x,y).")
      }
    }
    
    return(seq(bounds[1], bounds[2], length.out = k))
  }
  
  u <- function(x, zk, Tk, h_Tk, hprime_Tk) {
    j <- findInterval(x, zk) + 1
    calc_u <- function(x){
      u_result <- h_Tk[j] + (x-Tk[j])*hprime_Tk[j]
      return(u_result)
    }
    return(calc_u(x))
  }
  #   
  # l <- function(x, Tk, h_tk, hprime_Tk) {
  #   n <- length(Tk)
  #   calc_l <- function(x) {
  #     l_result <- ((Tk[2:n] - x) * h_Tk[1:(n-1)] + (x - Tk[1:(n-1)]) * h_Tk[2:n]) / (Tk[2:n] - Tk[1:(n-1)])
  #     return(l_result)
  #   }
  #   return(calc_l)
  # }
  
  
  l <- function(x, Tk, h_Tk, hprime_Tk) {
    j <- findInterval(x, Tk)
    if( (j == 0) || (j == length(Tk)) ) {
      return(-Inf)
    } else {
      calc_l <- function(x) {
        l_result <- ((Tk[j+1] - x) * h_Tk[j] + (x - Tk[j]) * h_Tk[j+1]) / (Tk[j+1] - Tk[j])
        return(l_result)
      }
      
      return(calc_l(x))
      
    }
  }
  
  ### outline: come back to this later
  s <- function(x, xk, Tk, h_Tk, hprime_Tk) {
    calc_s <- function(x){
      return(u_out / u_out_prime)
    }
    return(calc_s)
  }
  
  s_new <- function(x, zk, Tk, h_Tk, hprime_Tk) {
    j <- findInterval(x, zk) + 1
    s_numer <- function(x) {
      exp(h_Tk[j] + (x-Tk[j])*hprime_Tk[j])
    }
    s_int <- integrate(s_numer, lower = zk[j-1], upper = zk[j])$value
    
    s_Tk <- function(x){
      (exp(h_Tk[j] + (x-Tk[j])*hprime_Tk[j])) / s_int
    }
    return(s_Tk)
  }
  
  calc_probs <- function(Tk, zk, h_Tk, hprime_Tk) {
    
    z_all <- c(bounds[1], zk, bounds[2])
    z_plus_one <- z_all[2:length(z_all)]
    num_bins <- length(z_all) - 1
    cdf_vals <- rep(0, num_bins)
    z_start <- z_all[1:num_bins]
    
    cdf_vals <- exp(h_Tk)/hprime_Tk * (exp(z_plus_one*hprime_Tk) - exp(z_start*hprime_Tk))
    
    normalizingConstant <- sum(cdf_vals)
    
    return(cdf_vals/normalizingConstant)
  }
  
  sample_sk <- function(Tk, zk, h_Tk, hprime_Tk) {
    
    z_all <- c(bounds[1], zk, bounds[2])
    probs <- calc_probs(Tk, zk, h_Tk, hprime_Tk)
    
    j <- sample(length(probs), size = 1, prob = probs)
    u <- runif(1)
    x_star <- (1/hprime_Tk[j]) * log(u * (exp(hprime_Tk[j]*z_all[j+1]) - exp(hprime_Tk[j]*z_all[j]))
                                     + exp(hprime_Tk[j]*z_all[j]))
    
    return(x_star)
  }
  
  
############### FUNCTION START ###########################################
  check_input(f, n, bounds, k)
  
  n <- as.integer(n)
  x_init <- as.integer(x_init)
  k <- as.integer(k)
  
  
  h <- function(x) {
    
    logf <- log(f(x))
    assertthat::assert_that(is.finite(logf), msg = sprintf("Invalid Bounds case logf: %s.", x))
    return(log(f(x)))
  }
  
  
  hprime <- function(x) {
    der <- numDeriv::grad(h, x)
    return(der)
  }
  
  # initializing sample vector
  samps <- rep(NA, n)
  
  # INITIALIZING STEP
  
  Tk <- initialize_abcissae(k, bounds, hprime)
  
  num_samps <- 1
  
  h_Tk <- sapply(Tk, h)
  hprime_Tk <- sapply(Tk, hprime)
  check_log_concave(hprime_Tk)
  
  zk <- calc_z(Tk, h_Tk, hprime_Tk)
  
  # for debugging
  #else_count <- 0
  
  ### SAMPLING STEP
  while(num_samps <= 1000) {
    
    # real xstar
    #xstar <- sample_sk(Tk, zk, h_Tk, hprime_Tk)
    
    # temporary xstar
    xstar <- sample(bounds[1]:bounds[2], 1)
    
    #### check that xstar in bounds[1], bounds[2]
    
    w <- runif(1)
    
    # squeezing test
    if(w <= exp(l(xstar, Tk, h_Tk, hprime_Tk) - u(xstar, zk, Tk, h_Tk, hprime_Tk))) {
      
      samps[num_samps] <- xstar
      num_samps <- num_samps + 1
    }
    
    # rejection test
    else if (w <= exp(h(xstar) - u(xstar, zk, Tk, h_Tk, hprime_Tk))) {
      samps[num_samps] <- xstar
      num_samps <- num_samps + 1
    }
    
    ### UPDATING STEP
    else {
      #else_count <- else_count + 1
      h_xstar <- h(xstar)
      hprime_xstar <- hprime(xstar)
      Tk <- sort(c(Tk, xstar))
      
      # not sure if this stuff is right
      h_Tk <- append(h_Tk, h_xstar, after = (which(Tk == xstar) - 1))
      hprime_Tk <- append(hprime_Tk, hprime_xstar, after = (which(Tk == xstar) - 1))
      
      # for debugging later
      #print(else_count)
      
      #### construct u_k+1(x), s, l and increment k: k <- k + 1
      
      
      ### at some point in this loop, also check log-concavity
    }
    
    
  }
  
  return(samps)
  
}




f <- function(x, mean = 55, sigma = 1) {
  
  1/sqrt(2*pi*sigma^2)*exp(-1/2*(x-mean)^2/(sigma^2))
}

test <- ars(f = f, n = 1000, bounds = c(45, 55))

hist(test)
