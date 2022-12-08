
ars <- function(f, n = 1000, bounds = c(-Inf, Inf), x_init = 1, k = 20, ...) {
  
  assertthat::assert_that(is.function(f), msg = "f must be a function")
  assertthat::assert_that(is.numeric(n), msg = "n must be an integer")
  assertthat::assert_that(is.numeric(x_init), msg = "x_init must be numeric")
  assertthat::assert_that((is.numeric(k) && (k >= 2)), msg = "k must be an integer greater than or equal to 2")
  
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
  
  assertthat::assert_that((x_init >= bounds[1]) && (x_init <= bounds[2]), msg = "x_init must be inside bounds")
  
  
  ################################################################################
  
  check_log_concave <- function(x){
    # ensure that h'(x) is decreasing monotonically
    n <- length(x)
    assertthat::assert_that(sum(x[2:n] - x[1:n-1] <= 1e-8) == n-1, msg = "Function is not log-concave")
  }
  
  ########################################################################
  
  # initialize_abcissae <- function(x_init, k, hprime, bounds) {
  #   
  #   if((bounds[1] == -Inf) && (bounds[2] == Inf)) {
  #     x <- x_init - 10
  #     y <- x_init + 10
  #     
  #     while((!is.finite(hprime(x))) && (!is.finite(hprime(y)))) {
  #       x <- x + 1
  #       y <- y - 1
  #       
  #       if (x >= y) {
  #         stop("Invalid Bounds case (-Inf, Inf).")
  #       }
  #       
  #     }
  #     
  #     bounds[1] <- x
  #     bounds[2] <- y
  #     bounds[1] <<- x
  #     bounds[2] <<- y
  #     
  #   }
  #   
  #   if((bounds[1] == -Inf) && (bounds[2] != Inf)) {
  #     
  #     x <- bounds[2] - 10
  #     
  #     while((hprime(x) <= 0) && (hprime(x) != Inf) && (x < bounds[2])) {
  #       x <- x + 1
  #     }
  #     
  #     if((hprime(x) == Inf) | (x >= bounds[2])) {
  #       stop("Invalid Bounds case (-Inf, x).")
  #     }
  #     
  #     bounds[1] <- x
  #     bounds[1] <<- x
  #     
  #   }
  #   
  #   if((bounds[1] != -Inf) && (bounds[2] == Inf)) {
  #     
  #     y <- bounds[1] + 10
  #     
  #     while((hprime(y) >= 0) && (hprime(y) != -Inf) && (y > bounds[1])) {
  #       y <- y - 1
  #     }
  #     
  #     if((hprime(y) == -Inf) | (y <= bounds[1])) {
  #       stop("Invalid Bounds case (x, Inf).")
  #     }
  # 
  #     bounds[2] <- y
  #     bounds[2] <<- y
  #     
  #   }
  #   
  #   if((bounds[1] != -Inf) && (bounds[2] != Inf)) {
  #     
  #     if((!is.finite(hprime(bounds[1]))) || (!is.finite(hprime(bounds[2])))) {
  #       stop("Function f not defined on bounds case (x,y).")
  #     }
  #   }
  #   
  #   return(seq(bounds[1], bounds[2], length.out = k))
  # }
  
  
  initialize_abcissae <- function(x_init, k, hprime, bounds) {
    x1 <- bounds[1]
    xk <- bounds[2]
    #print(paste("xinit:", x_init))
    #print(paste("Boundsss 1:", bounds[1]))
    #print(paste("Bounds2:", bounds[2]))
    inc <- 1
    if (bounds[1] == -Inf) {
      x1 <- x_init
      while (hprime(x1) <= 0){
        x1 <- x1 - inc
        inc <- inc * 2
      }
      
    }
    if (bounds[2] == Inf) {
      x2 <- x_init
      while (hprime(xk) >= 0) {
        xk <- xk + inc
        inc <- inc * 2
      }
      
    }
    
    if ((bounds[1] != -Inf) && (bounds[2] != Inf)) {
      print("Finite bounds")
      x1 <- bounds[1]
      xk <- bounds[2]
      #print(paste("Bound 1:", x1))
      #print(paste("Bound 2:", x2))
    }
    #print(paste("x1:", x1))
    #print(paste("xk:", xk))
    return(seq(x1, xk, length.out = k))
  }
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  calc_z <- function(tk, h_tk, hprime_tk) {
    n <- length(tk)
    if (var(hprime_tk) <= 1e-8) {
      return(rep(0, n-1))
    }
    
    return((h_tk[2:n] - h_tk[1:(n-1)] - tk[2:n] * hprime_tk[2:n] + tk[1:(n-1)]* hprime_tk[1:(n-1)])/(hprime_tk[1:(n-1)]-hprime_tk[2:n]))
    
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
  
  # calc_probs <- function(Tk, zk, h_Tk, hprime_Tk) {
  #   print(zk)
  #   
  #   z_all <- c(bounds[1], zk, bounds[2])
  #   num_bins <- length(z_all) - 1
  #   
  #   z_2 <- z_all[2:length(z_all)]
  #   z_1 <- z_all[1:num_bins]
  #   
  #   s <- exp(h_Tk - hprime_Tk*Tk)
  #   s[is.infinite(s) | is.na(s)] <- 0
  #   
  #   cdf_vals <- c()
  #   cdf_vals[s!=0] <- (s/hprime_Tk) * (exp(z_2*(hprime_Tk)) - exp(z_1*(hprime_Tk))) 
  #   cdf_vals[s==0] <- s[s==0]*(z_2[s==0] - z_1)
  #   
  #   normalizing_constant <- sum(cdf_vals)
  #   normalized_b <- s/normalizing_constant
  #   
  #   weighted_prob <- c()
  #   weighted_prob[s!=0] <- normalized_b/hprime_Tk*(exp(hprime_Tk*z_2) - exp(hprime_Tk*z_1))
  #   weighted_prob[s==0] <- normalized_b*(z_2-z_1)
  #   return(list(normalized_b, weighted_prob))
  # }
  
  
  # calc_probs <- function(Tk, zk, h_Tk, hprime_Tk) {
  #   
  #   
  #   z_all <- c(bounds[1], zk, bounds[2])
  #   num_bins <- length(z_all) - 1
  #   cdf_vals <- rep(0, num_bins)
  #   
  #   z_plus_one <- z_all[2:length(z_all)]
  #   z_start <- z_all[1:num_bins]
  #   
  #   s <- exp(h_Tk - hprime_Tk*Tk)
  #   s[which(is.infinite(s))] <- 0
  #   s[which(is.na(s))] <- 0
  #   
  #   cdf_vals <- ifelse(s==0, (s/hprime_Tk) * (exp(z_plus_one*(hprime_Tk)) - exp(z_start*(hprime_Tk))), 
  #                      z_plus_one-z_start)
  #   
  #   normalizing_constant <- sum(cdf_vals)
  #   normalized_b <- s/normalizing_constant
  #   
  #   weighted_prob <- c()
  #   weighted_prob[s!=0] <- normalized_b/hprime_Tk*(exp(hprime_Tk*z_plus_one) - exp(hprime_Tk*z_start))
  #   weighted_prob[s==0] <- normalized_b*(z_plus_one-z_start)
  #   return(list(normalized_b, weighted_prob))
  # }
  
  calc_probs <- function(Tk, zk, h_Tk, hprime_Tk) {
    
    
    z_all <- c(bounds[1], zk, bounds[2])
    num_bins <- length(z_all) - 1
    
    z_2 <- z_all[2:length(z_all)]
    z_1 <- z_all[1:num_bins]
    
    s <- exp(h_Tk - hprime_Tk*Tk)
    s[is.infinite(s) | is.na(s)] <- 0
    
    cdf_vals <- s*(z_2 - z_1)
    cdf_vals[hprime_Tk!=0] <- (s[hprime_Tk!=0]/hprime_Tk[hprime_Tk!=0]) * (exp(z_2[hprime_Tk!=0]*(hprime_Tk[hprime_Tk!=0])) - exp(z_1[hprime_Tk!=0]*(hprime_Tk[hprime_Tk!=0]))) 
    #cdf_vals[hprime_Tk==0] <- s*(z_2 - z_1)
    
    normalizing_constant <- sum(cdf_vals)
    normalized_b <- s/normalizing_constant
    
    weighted_prob <- normalized_b*(z_2-z_1)
    weighted_prob[hprime_Tk!=0] <- normalized_b[hprime_Tk!=0]/hprime_Tk[hprime_Tk!=0]*(exp(hprime_Tk[hprime_Tk!=0]*z_2[hprime_Tk!=0]) - exp(hprime_Tk[hprime_Tk!=0]*z_1[hprime_Tk!=0]))
    weighted_prob[is.infinite(weighted_prob) | weighted_prob <= 0] <- 0
    return(list(normalized_b, weighted_prob))
  }
  
  
  
  
  
  sample_sk <- function(Tk, zk, h_Tk, hprime_Tk) {
    
    z_all <- c(bounds[1], zk, bounds[2])
    
    b <- calc_probs(Tk, zk, h_Tk, hprime_Tk)[[1]]
    prob <- calc_probs(Tk, zk, h_Tk, hprime_Tk)[[2]]
    prob[(prob <= 0) | (is.na(prob))] <- 0
    #print(b)
    #print(prob)
    i <- sample(length(prob), size = 1, prob = prob)
    u <- runif(1)
    x_star <- ifelse(hprime_Tk[i] == 0, z_all[i]+prob[i]*u/b[i],
                     (1/hprime_Tk[i]) * log(u * (hprime_Tk[i]*prob[i]/b[i])+ exp(hprime_Tk[i]*z_all[i])))
    
    return(x_star)
  }
  
  
  
############### FUNCTION START ###########################################
  
  n <- as.integer(n)
  k <- as.integer(k)
  
  
  h <- function(x) {
    
    logf <- log(f(x, ...))
    #assertthat::assert_that(is.finite(logf), msg = sprintf("Invalid Bounds case logf: %s.", x))
    return(logf)
  }
  
  
  hprime <- function(x) {
    der <- numDeriv::grad(h, x)
    #der <- (h(x + 1e-8) - h(x))/1e-8
    return(der)
  }
  
  # initializing sample vector
  samps <- rep(NA, n)
  
  # INITIALIZING STEP
  
  Tk <- initialize_abcissae(x_init, k, hprime, bounds)
  print(Tk)
  
  num_samps <- 1
  print("Found absissae")
  h_Tk <- sapply(Tk, h)
  print("Found h")
  print(h_Tk)
  hprime_Tk <- sapply(Tk, hprime)
  print("Found derivative")
  print(hprime_Tk)
  hprime_Tk[is.na(hprime_Tk)]
  hprime_Tk[is.na(hprime_Tk)] <- 0
  print(hprime_Tk)
  print(var(hprime_Tk))
  check_log_concave(hprime_Tk)
  
  zk <- calc_z(Tk, h_Tk, hprime_Tk)
  zk[is.na(zk)] <- 0
  print(zk)
  
  
  ### SAMPLING STEP
  while(num_samps <= n) {
    
    # sample xstar from sk
    xstar <- sample_sk(Tk, zk, h_Tk, hprime_Tk)
    print(xstar)
    print(paste("l:", l(xstar, Tk, h_Tk, hprime_Tk) ))
    print(paste("h:", h(xstar)))
    print(zk)
    print(Tk)
    print(paste("u", u(xstar, zk, Tk, h_Tk, hprime_Tk)))
    assertthat::assert_that((xstar >= bounds[1]) && (xstar <= bounds[2]), msg = "xstar not in bounds")
    #zk <- sort(zk)
    assertthat::assert_that((l(xstar, Tk, h_Tk, hprime_Tk) <= h(xstar)) && (h(xstar) <= u(xstar, zk, Tk, h_Tk, hprime_Tk)), msg = "lhu test: Not log concave")
    w <- runif(1)
    
    # squeezing test
    if(w <= exp(l(xstar, Tk, h_Tk, hprime_Tk) - u(xstar, zk, Tk, h_Tk, hprime_Tk))) {
      
      samps[num_samps] <- xstar
      num_samps <- num_samps + 1
    }
    
    # rejection test
    else {
      print("rejected")
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

hist(ars(f = dgamma, n = 1000, bounds = c(1,10), shape = 9), freq = F)
#ars(f = dunif, n = 1000, bounds = c(0,1))
# curve(dnorm(x, 0, 1), 0, 1, add = TRUE, col = "red")
# 
# test <- ars(f = dnorm, n = 1000, bounds = c(55, 60),  x_init = 47, mean = 50)
# hist(test, freq = F)
# curve(dnorm(x, 50, 1), 55, 60,  add = TRUE, col = "red")
# 
# 
# # testing with gamma
#test <- ars(f = dgamma, n = 1000,  bounds = c(0,10), x_init = 5, shape = 9, rate = 2)
# 
# hist(test, freq = F)
# curve(dgamma(x, 9, 2), 1, 10, add = TRUE, col = "red")
# 
# 
# # testing with logistics
# test <- ars(f = dlogis, n = 1000,  bounds = c(1,5), x_init = 0.5)
# hist(test, freq = F)
# curve(dlogis(x, 0, 1), 0, 10, add = TRUE, col = "red")
# 
# 
# # testing with unif
# test <- ars(f = dnorm, n = 1000,  bounds = c(0,Inf), x_init = 1)
# 
# test <- ars(f = dnorm, n = 1000, bounds = c(10, 15), min = 10, max = 15, x_init = 11)
# hist(test, freq = F)

# test <- ars(f = dnorm, n = 1000, bounds = c(0,1), x_init = 0.5)
# hist(test, freq = F)


test <- ars(f = dunif, n = 1000, bounds = c(10,15), x_init = 11, min = 10, max = 15)
hist(test, freq = F)
curve(dbeta(x, 3, 4), 0.01, .99, add = TRUE, col = "red")




####
