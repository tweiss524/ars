
##### AUXILIARY FUNCTIONS FOR ARS #####
#######################################


# function check_log_concave takes vector of derivatives as input
# and checks that it is monotonically decreasing
check_log_concave <- function(x) {
  n <- length(x)
  assertthat::assert_that(sum(x[2:n] - x[1:(n-1)] <= 1e-8) == n-1, msg = "Function is not log-concave")
}

# function initialize_abscissae take in x_init, hprime function, and bounds as input
# uses x_init as a reference within the bounds to produce a vector of length
# 20, representing the initial abscissae points
initialize_abscissae <- function(x_init, hprime, bounds) {

  x1 <- bounds[1]
  xk <- bounds[2]
  inc <- 0.25

  # case where left bound is infinite
  if (x1 == -Inf) {
    x1 <- x_init - 0.1

    # run until h'(x) > 0, each time exponentially increasing the step size
    while(!is.na(hprime(x1)) && hprime(x1) <= 0){
      x1 <- x1 - inc
      inc <- inc * 2
    }

    if(is.na(hprime(x1))) {
      stop("Please provide x_init with a finite derivative")
    }
  }

  # case where right bound is infinite
  if (xk == Inf) {
    xk <- x_init + 0.1

    # run until h'(x) < 0, each time exponentially increasing the step size
    while (!is.na(hprime(xk)) && hprime(xk) >= 0) {
      xk <- xk + inc
      inc <- inc * 2
    }

    if(is.na(hprime(xk))) {
      stop("Please provide x_init with a finite derivative")
    }
  }

  # case where both bounds are finite
  if ((bounds[1] != -Inf) && (bounds[2] != Inf)) {
    x1 <- bounds[1]
    xk <- bounds[2]
  }
  return(seq(x1, xk, length.out = 20))
}


# function calc_z takes in Tk, h_Tk, and hprime_Tk to calculate zk,
# the intersection points of the tangents for each elements in Tk.
calc_z <- function(tk, h_tk, hprime_tk) {
  n <- length(tk)
  
  # handles single value constant derivative case
  if (n == 1) {
    return(c())
  }
  
  z_res <- (h_tk[2:n] - h_tk[1:(n-1)] - tk[2:n] * hprime_tk[2:n] + tk[1:(n-1)] *
              hprime_tk[1:(n-1)])/(hprime_tk[1:(n-1)] - hprime_tk[2:n])
  
  return(z_res)

}


# u function takes in a sampled x* value, zk (interection points).
# Tk (abscissae points), h(Tk), and h'(Tk), and 
# calculates the rejection envelope on Tk
u <- function(x, zk, Tk, h_Tk, hprime_Tk) {
  
  # handle the case where zk == NULL, meaning no intersection points
  if (length(zk) == 0){
    j <- 1
  } else {
    j <- findInterval(x, zk) + 1
  }
  
  calc_u <- function(x) {
    u_result <- h_Tk[j] + (x - Tk[j]) * hprime_Tk[j]
    return(u_result)
  }
  
  return(calc_u(x))
}


# l function takes in a sampled x* value, zk (interection points).
# Tk (abscissae points), h(Tk), and h'(Tk), and 
# calculates the squeezing function on Tk
l <- function(x, zk, Tk, h_Tk, hprime_Tk) {
  
  # handle the case where there is only one initial point
  if (length(Tk) == 1) {
    return(u(x, zk, Tk, h_Tk, hprime_Tk))
  } else {
    j <- findInterval(x, Tk)
    
    # j < first element and > last element
    if( (j == 0) || (j == length(Tk)) ) {
      return(-Inf)
    } else {
      
      calc_l <- function(x) {
        l_result <- ((Tk[j+1] - x) * h_Tk[j] + (x - Tk[j]) * h_Tk[j+1]) /
          (Tk[j+1] - Tk[j])
        
        return(l_result)
        
      }
      return(calc_l(x))
    }
  }
}

# function calc_probs takes as input initial abscissae points Tk, 
# intersection points zk, h(Tk), h'(Tk), and user-defined bounds
# and returns a probability vector 
calc_probs <- function(Tk, zk, h_Tk, hprime_Tk, bounds) {

  # appending bounds to zk
  z_all <- c(bounds[1], zk, bounds[2])
  num_bins <- length(z_all) - 1

  z_2 <- z_all[2:length(z_all)]
  z_1 <- z_all[1:num_bins]
  z_ind <- which(hprime_Tk != 0)

  u_z1 <- u(z_1, zk, Tk, h_Tk, hprime_Tk)
  u_z2 <- u(z_2, zk, Tk, h_Tk, hprime_Tk)

  unnormalized_prob <- exp(u_z1) * (z_2 - z_1)
  unnormalized_prob[z_ind] <- (exp(u_z2[z_ind]) - exp(u_z1[z_ind])) / hprime_Tk[z_ind]

  # normalize the probability
  normalized_prob <- unnormalized_prob/sum(unnormalized_prob)



  return(list(normalized_prob,unnormalized_prob))
}

# function sample_sk takes as input initial abscissae points Tk, 
# intersection points zk, h(Tk), h'(Tk), and user-defined bounds, and
# computes the inverse CDF for sampling, giving probability weights
# based on the probability vector returned from calc_probs()
sample_sk <- function(Tk, zk, h_Tk, hprime_Tk, bounds) {

  z_all <- c(bounds[1], zk, bounds[2])
  num_bins <- length(z_all) - 1
  z_1 <- z_all[1:num_bins]

  p <- calc_probs(Tk, zk, h_Tk, hprime_Tk, bounds)
  
  # normalized
  prob <- p[[1]]
  
  # unnormalized
  unnorm_prob <- p[[2]]
  #prob[(prob <= 0) | (is.na(prob))] <- 0
  u_z1 <- u(z_1, zk, Tk, h_Tk, hprime_Tk)
  i <- sample(length(prob), size = 1, prob = prob)
  unif <- runif(1)
  
  # sampled point from function
  x_star <- ifelse(hprime_Tk[i] == 0, z_all[i] + unif*(z_all[i+1]-z_all[i]),
                   (log(unif * unnorm_prob[i] * hprime_Tk[i] + exp(u_z1[i])) -
                      h_Tk[i] + hprime_Tk[i] * Tk[i]) / hprime_Tk[i])
  return(x_star)
}



