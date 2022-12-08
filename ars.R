
ars <- function(f, n = 1000,
                bounds = c(-Inf, Inf), k = 20,
                x_init = 1, ...) {
  
  assertthat::assert_that(is.function(f), msg = "f must be a function")
  assertthat::assert_that(is.numeric(n), msg = "n must be an integer")
  assertthat::assert_that(is.numeric(x_init), msg = "x_init must be an integer")
  
  assertthat::assert_that(is.vector(bounds) && (length(bounds) == 2) && (is.numeric(bounds)), msg = "Bounds must be numeric vector of length 2")
  
  if (bounds[1] > bounds[2]) {
    
    bounds <- c(bounds[2], bounds[1])
    
    warning(paste0("Lower bound must be smaller than upper bound. ", 
                   sprintf("New bounds are (%s, %s)", bounds[1], bounds[2])))
    
  }
  
  if (bounds[1] == bounds[2]) {
    
    bounds <- c(-Inf, Inf)
    
    warning(paste0("Lower bound must be smaller than upper bound. ", 
                   sprintf("New bounds are (%s, %s)", bounds[1], bounds[2])))
    
  }
  
  assertthat::assert_that((is.numeric(k)) && (k > 0), msg = "k must be a positive integer")
  
  assertthat::assert_that((x_init >= bounds[1]) && (x_init <= bounds[2]), msg = "x_init must be inside bounds")
  
  
  ################################################################################
  
  check_log_concave <- function(x){
    ep <- 1e-8
    n <- length(x)
    assertthat::assert_that(sum(x[2:n] - x[1:n-1] < ep) == n-1, msg = "Function is not log-concave")
  }
  
  
  ###############################################################################
  
  calc_z <- function(tk, h_tk, hprime_tk) {
    n <- length(tk)
    
    return((h_tk[2:n] - h_tk[1:(n-1)] - tk[2:n] * hprime_tk[2:n] + tk[1:(n-1)]* hprime_tk[1:(n-1)])/(hprime_tk[1:(n-1)]-hprime_tk[2:n]))
    
  }
  
  ########################################################################
  
  initialize_abcissae <- function(x_init, k, hprime, bounds) {
    
    if((bounds[1] == -Inf) && (bounds[2] == Inf)) {
      x <- x_init - 10
      y <- x_init + 10
      
      while((!is.finite(hprime(x))) && (!is.finite(hprime(y)))) {
        x <- x + 1
        y <- y - 1
        
        if (x >= y) {
          stop("Invalid Bounds case (-Inf, Inf).")
        }
        
      }
      
      bounds[1] <- x
      bounds[2] <- y
      bounds[1] <<- x
      bounds[2] <<- y
      
    }
    
    if((bounds[1] == -Inf) && (bounds[2] != Inf)) {
      
      x <- bounds[2] - 10
      
      while((hprime(x) <= 0) && (hprime(x) != Inf) && (x < bounds[2])) {
        x <- x + 1
      }
      
      if((hprime(x) == Inf) | (x >= bounds[2])) {
        stop("Invalid Bounds case (-Inf, x).")
      }
      
      bounds[1] <- x
      bounds[1] <<- x
      
    }
    
    if((bounds[1] != -Inf) && (bounds[2] == Inf)) {
      
      y <- bounds[1] + 10
      
      while((hprime(y) >= 0) && (hprime(y) != -Inf) && (y > bounds[1])) {
        y <- y - 1
      }
      
      if((hprime(y) == -Inf) | (y <= bounds[1])) {
        stop("Invalid Bounds case (x, Inf).")
      }
  
      bounds[2] <- y
      bounds[2] <<- y
      
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
  
  calc_probs <- function(Tk, zk, h_Tk, hprime_Tk) {
    
    
    z_all <- c(bounds[1], zk, bounds[2])
    num_bins <- length(z_all) - 1
    
    z_2 <- z_all[2:length(z_all)]
    z_1 <- z_all[1:num_bins]
    
    s <- exp(h_Tk - hprime_Tk*Tk)
    s[is.infinite(s) | is.na(s)] <- 0
    
    cdf_vals <- c()
    cdf_vals[s!=0] <- (s/hprime_Tk) * (exp(z_2*(hprime_Tk)) - exp(z_1*(hprime_Tk))) 
    cdf_vals[s!=0] <- z_2 - z_1
    
    normalizing_constant <- sum(cdf_vals)
    normalized_b <- s/normalizing_constant
    
    weighted_prob <- c()
    weighted_prob[s!=0] <- normalized_b/hprime_Tk*(exp(hprime_Tk*z_2) - exp(hprime_Tk*z_1))
    weighted_prob[s==0] <- normalized_b*(z_2-z_1)
    return(list(normalized_b, weighted_prob))
  }
  
  sample_sk <- function(Tk, zk, h_Tk, hprime_Tk) {
    
    z_all <- c(bounds[1], zk, bounds[2])
    
    b <- calc_probs(Tk, zk, h_Tk, hprime_Tk)[[1]]
    prob <- calc_probs(Tk, zk, h_Tk, hprime_Tk)[[2]]
    prob[(prob <= 0) | (is.na(prob))] <- 0
    
    i <- sample(length(prob), size = 1, prob = prob)
    u <- runif(1)
    x_star <- ifelse(hprime_Tk[i] == 0, z_all[i]+prob[i]*u/b[i],
                     (1/hprime_Tk[i]) * log(u * (hprime_Tk[i]*prob[i]/b[i])+ exp(hprime_Tk[i]*z_all[i])))
    
    return(x_star)
  }
  
  
  
############### FUNCTION START ###########################################
  
  n <- as.integer(n)
  x_init <- as.integer(x_init)
  k <- as.integer(k)
  
  
  h <- function(x) {
    
    logf <- log(f(x, ...))
    assertthat::assert_that(is.finite(logf), msg = sprintf("Invalid Bounds case logf: %s.", x))
    return(logf)
  }
  
  
  hprime <- function(x) {
    der <- numDeriv::grad(h, x)
    return(der)
  }
  
  # initializing sample vector
  samps <- rep(NA, n)
  
  # INITIALIZING STEP
  
  Tk <- initialize_abcissae(x_init, k, hprime, bounds)
  
  num_samps <- 1
  
  h_Tk <- sapply(Tk, h)
  hprime_Tk <- sapply(Tk, hprime)
  check_log_concave(hprime_Tk)
  
  zk <- calc_z(Tk, h_Tk, hprime_Tk)
  
  
  ### SAMPLING STEP
  while(num_samps <= n) {
    
    # sample xstar from sk
    xstar <- sample_sk(Tk, zk, h_Tk, hprime_Tk)
    
    assertthat::assert_that((xstar >= bounds[1]) && (xstar <= bounds[2]), msg = "xstar not in bounds")
    zk <- sort(zk)
    assertthat::assert_that((l(xstar, Tk, h_Tk, hprime_Tk) <= h(xstar)) && (h(xstar) <= u(xstar, zk, Tk, h_Tk, hprime_Tk)), msg = "lhu test: Not log concave")
    w <- runif(1)
    
    # squeezing test
    if(w <= exp(l(xstar, Tk, h_Tk, hprime_Tk) - u(xstar, zk, Tk, h_Tk, hprime_Tk))) {
      
      samps[num_samps] <- xstar
      num_samps <- num_samps + 1
    }
    
    # rejection test
    else {
      
      h_xstar <- h(xstar)
      hprime_xstar <- hprime(xstar)
      
      if (w <= exp(h_xstar - u(xstar, zk, Tk, h_Tk, hprime_Tk))) {
        samps[num_samps] <- xstar
        num_samps <- num_samps + 1
      }
      
      Tk <- sort(c(Tk, xstar))
      
      h_Tk <- append(h_Tk, h_xstar, after = (which(Tk == xstar) - 1))
      hprime_Tk <- append(hprime_Tk, hprime_xstar, after = (which(Tk == xstar) - 1))
      
      check_log_concave(hprime_Tk)
      
      zk <- calc_z(Tk, h_Tk, hprime_Tk)
      
      } # end else
    } # end while
    
  return(samps)
}


# testing with normal
test <- ars(f = dnorm, n = 1000, bounds = c(45,55),  x_init = 47, mean = 50)

hist(test, freq = F)
curve(dnorm(x, 50, 1), 45, 55, add = TRUE, col = "red")


# testing with gamma
test <- ars(f = dgamma, n = 1000,  bounds = c(1,10), x_init = 5, shape = 9, rate = 2)

hist(test, freq = F)
curve(dgamma(x, 9, 2), 1, 10, add = TRUE, col = "red")


# testing with logistics
test <- ars(f = dlogis, n = 1000,  bounds = c(0,10), x_init = 0.5)

# testing with unif
test <- ars(f = dunif, n = 1000,  bounds = c(0,1), x_init = 0.5)


