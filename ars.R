
ars <- function(f, n = 1000, bounds = c(-Inf, Inf), x_init = 1, k = 20, ...) {
  #op <- optim(1, f, method = 'BFGS', control = list(fnscale=-1))
  #x_init
  #print("Val:")
  #print(op$value)
  #print("X_init:")
  #print(op$par)
  #print("f(xinit)")
  #print(f(op$par))
  #if (is.na(x_init)) {x_init <- 1}
  assertthat::assert_that(is.function(f), msg = "f must be a function")
  assertthat::assert_that(is.numeric(n), msg = "n must be an integer")
  assertthat::assert_that(is.numeric(x_init), msg = "x_init must be numeric")
  assertthat::assert_that((is.numeric(k) && (k >= 3)), msg = "k must be an integer greater than or equal to 3")
  
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
  
  assertthat::assert_that((x_init >= bounds[1]) && (x_init <= bounds[2]), msg = "x_init points must be inside bounds")
  
  ################################################################################
  
  check_log_concave <- function(x){
    # ensure that h'(x) is decreasing monotonically
    n <- length(x)
    assertthat::assert_that(sum(x[2:n] - x[1:n-1] <= 1e-8) == n-1, msg = "Function is not log-concave")
  }
  
  ########################################################################
  
#   
#   initialize_abcissae <- function(x_init, k, hprime, bounds) {
#     Tk <- c(x_init)
#     
#     x1 <- bounds[1]
#     xk <- bounds[2]
#     inc <- 0.25
#     hpTk <- sapply(Tk, hprime)
#     if (bounds[1] == -Inf) {
#       assertthat::assert_that(!all(hpTk < 0), msg = "h'(x_init) must be > 0")
#       x1 <- which(hpTk > 0)[1]
#       }
# 
#     if (bounds[2] == Inf) {
#       assertthat::assert_that(!all(hpTk > 0), msg = "h'(x_init) must be < 0")
#       xk <- sort(which(hpTk < 0), decreasing = T)[1]
#       }
#     
#     if ((bounds[1] != -Inf) && (bounds[2] != Inf)) {
#       print("Finite bounds")
#       x1 <- bounds[1]
#       xk <- bounds[2]
#       #print(paste("Bound 1:", x1))
#       #print(paste("Bound 2:", x2))
#     }
#     #print(paste("x1:", x1))
#     #print(paste("xk:", xk))
#     return(seq(x1, xk, length.out = 2))
# }

  
  
  initialize_abcissae <- function(x_init, k, hprime, bounds) {
    x1 <- bounds[1]
    xk <- bounds[2]
    #print(paste("xinit:", x_init))
    #print(paste("Boundsss 1:", bounds[1]))
    #print(paste("Bounds2:", bounds[2]))
    inc <- 0.25
    if (bounds[1] == -Inf) {
      x1 <- x_init
      while (hprime(x1) <= 0){
        x1 <- x1 - inc
        inc <- inc * 2
      }
      
    }
    if (bounds[2] == Inf) {
      xk <- x_init
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
    if (n == 1) {
      return(c())
    }
    
    return((h_tk[2:n] - h_tk[1:(n-1)] - tk[2:n] * hprime_tk[2:n] + tk[1:(n-1)]* hprime_tk[1:(n-1)])/(hprime_tk[1:(n-1)]-hprime_tk[2:n]))
    
  }
  
  
  u <- function(x, zk, Tk, h_Tk, hprime_Tk) {
    if (length(zk)==0){
      j <- 1
    }else{
      j <- findInterval(x, zk) + 1
    }
    calc_u <- function(x){
      u_result <- h_Tk[j] + (x-Tk[j])*hprime_Tk[j]
      return(u_result)
    }
    return(calc_u(x))
  }
  
  
  l <- function(x, Tk, h_Tk, hprime_Tk) {
    if(length(Tk)==1){
      return(u(x, zk, Tk, h_Tk, hprime_Tk))
    }else{
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
  }
  
  calc_probs <- function(Tk, zk, h_Tk, hprime_Tk) {


    z_all <- c(bounds[1], zk, bounds[2])
    num_bins <- length(z_all) - 1
    print(paste("z-all:", z_all))

    z_2 <- z_all[2:length(z_all)]
    z_1 <- z_all[1:num_bins]

    s <- exp(h_Tk - hprime_Tk*Tk)
    s[is.infinite(s) | is.na(s)] <- 0
    print(paste("s", s))

    cdf_vals <- s*(z_2 - z_1)
    cdf_vals[hprime_Tk!=0] <- (s[hprime_Tk!=0]/hprime_Tk[hprime_Tk!=0]) * (exp(z_2[hprime_Tk!=0]*(hprime_Tk[hprime_Tk!=0])) - exp(z_1[hprime_Tk!=0]*(hprime_Tk[hprime_Tk!=0])))
    #cdf_vals[hprime_Tk==0] <- s*(z_2 - z_1)
    print(paste("cdf_vals:", cdf_vals))

    normalizing_constant <- sum(cdf_vals)
    normalized_b <- s/normalizing_constant
    print(paste("normalized_b:", normalized_b))

    weighted_prob <- normalized_b*(z_2-z_1)
    print(paste("weighted_prob:", weighted_prob))
    weighted_prob[hprime_Tk!=0] <- normalized_b[hprime_Tk!=0]/hprime_Tk[hprime_Tk!=0]*(exp(hprime_Tk[hprime_Tk!=0]*z_2[hprime_Tk!=0]) - exp(hprime_Tk[hprime_Tk!=0]*z_1[hprime_Tk!=0]))
    weighted_prob[is.infinite(weighted_prob) | weighted_prob <= 0] <- 0
    return(list(normalized_b, weighted_prob))
  }
  


  sample_sk <- function(Tk, zk, h_Tk, hprime_Tk) {

    z_all <- c(bounds[1], zk, bounds[2])

    b <- calc_probs(Tk, zk, h_Tk, hprime_Tk)[[1]]
    prob <- calc_probs(Tk, zk, h_Tk, hprime_Tk)[[2]]
    #prob[(prob <= 0) | (is.na(prob))] <- 0
    #print(b)
    #print("Probs:")
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
  print("Tk:")
  print(Tk)
  
  num_samps <- 1
  print("Found absissae")
  h_Tk <- sapply(Tk, h)
  print("Found h")
  print(h_Tk)
  
  # test if defined inside bounds by removing infinite log(f(x))
  Tk <- Tk[is.finite(h_Tk)]
  h_Tk <- h_Tk[is.finite(h_Tk)]
  assertthat::assert_that(length(h_Tk) > 0, msg = "Function not defined in bounds")
  
  
  hprime_Tk <- sapply(Tk, hprime)
  print("Found derivative")
  print(hprime_Tk)
  #print(hprime_Tk)
  #hprime_Tk[is.na(hprime_Tk)] <- 0
  h_Tk <- h_Tk[!is.na(hprime_Tk)]
  Tk <- Tk[!is.na(hprime_Tk)]
  hprime_Tk <- hprime_Tk[!is.na(hprime_Tk)]
  # 
  #assertthat::assert_that(length(h_Tk) > 0, msg = "Function not defined in bounds")
  print("hprime_Tk:")
  print(hprime_Tk)
  #print(sum(abs(hprime_Tk[2:length(hprime_Tk)] - hprime_Tk[1:(length(hprime_Tk) -1)]) <=  1e-8)  == (length(hprime_Tk)-1))
  
  # check if all h'(Tk) are the same and if they are, only keep one element
  len_hptk <- length(hprime_Tk)
  if (sum(abs(hprime_Tk[2:len_hptk] - hprime_Tk[1:(len_hptk-1)]) <=  1e-8)  == (len_hptk-1)) {
    print("if passed")
    hprime_Tk <- hprime_Tk[2]
    Tk <- Tk[2]
    h_Tk <- h_Tk[2]
  }
  # else if (length(unique(hprime_Tk)) == 2 && c(1,-1) %in% hprime_Tk) {
  #   hprime_Tk <- unique(hprime_Tk)
  #   Tk <- Tk[which(hprime_Tk == unique(hprime_Tk))]
  #   h_Tk <- h_Tk[which(hprime_Tk == unique(hprime_Tk))]
  # }

  check_log_concave(hprime_Tk)
  
  
  
  
  
  # Tk <- Tk[is.finite(h_Tk)]
  # h_Tk <- h_Tk[is.finite(h_Tk)]
  # hprime_Tk <- hprime_Tk[is.finite(h_Tk)]
  
  zk <- calc_z(Tk, h_Tk, hprime_Tk)
  #zk[is.na(zk)] <- 0
  print("zk:")
  print(zk)
  
  
  ### SAMPLING STEP
  while(num_samps <= n) {
    
    # sample xstar from sk
    xstar <- sample_sk(Tk, zk, h_Tk, hprime_Tk)
    #print(paste("xstar:", xstar))
    #print(paste("l:", l(xstar, Tk, h_Tk, hprime_Tk) ))
    #print(paste("h:", h(xstar)))
    # print(zk)
    # print(Tk)
    #print(paste("u", u(xstar, zk, Tk, h_Tk, hprime_Tk)))
    assertthat::assert_that((xstar >= bounds[1]) && (xstar <= bounds[2]), msg = "sampled xstar not in bounds")
    #zk <- sort(zk)
    #assertthat::assert_that((l(xstar, Tk, h_Tk, hprime_Tk) <= h(xstar)) && (h(xstar) <= u(xstar, zk, Tk, h_Tk, hprime_Tk)), msg = "lhu test: Not log concave")
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
      print("hstar:")
      print(h_xstar)
      if(!is.finite(h_xstar)) { break}
      hprime_xstar <- hprime(xstar)
      
      
      
      if (w <= exp(h_xstar - u(xstar, zk, Tk, h_Tk, hprime_Tk))) {
        samps[num_samps] <- xstar
        num_samps <- num_samps + 1
        print("accepted")
      }
      
      Tk <- sort(c(Tk, xstar))
      
      h_Tk <- append(h_Tk, h_xstar, after = (which(Tk == xstar) - 1))
      hprime_Tk <- append(hprime_Tk, hprime_xstar, after = (which(Tk == xstar) - 1))
      
      h_Tk <- h_Tk[!is.na(hprime_Tk)]
      Tk <- Tk[!is.na(hprime_Tk)]
      hprime_Tk <- hprime_Tk[!is.na(hprime_Tk)]
      
      
      
      check_log_concave(hprime_Tk)
      
      zk <- calc_z(Tk, h_Tk, hprime_Tk)
      print("zk")
      print(zk)
      
      } # end else
    } # end while
    
  return(samps)
}


# # testing with normal
# 
# hist(ars(f = dgamma, n = 1000, k = 3,  x_init = -1, bounds = c(-Inf, 0), shape = 3, rate = 4), freq = F)
# curve(dgamma(x, 3, 4), 0, add = TRUE, col = "red")
# ars(f = dunif, n = 1000, bounds = c(0,1))
# 
# 
#test <- ars(f = dnorm, n = 1000, mean = 10000, x_init = 10000)

#hist(ars(dcauchy, 1000, x_init = 0, bounds = c(-5,5)))
# hist(test, freq = F)
# curve(dnorm(x, 50, 1), 40, 60,  add = TRUE, col = "red")
# 
# # 
# # # testing with gamma
# test <- ars(f = dgamma, n = 1000,  bounds = c(1,10), x_init = 5, shape = 9, rate = 2)
# #
# hist(test, freq = F)
# curve(dgamma(x, 9, 2), 1, 10, add = TRUE, col = "red")
# 
# # 
# # # testing with logistics
#test <- ars(f = dlogis, n = 1000,  bounds = c(1,5), x_init = 2)
#hist(test, freq = F)
#curve(dlogis(x, 0, 1), 1, 5, add = TRUE, col = "red")
# 
# # 
# test <- ars(f = dnorm, n = 1000,  bounds = c(-Inf,0), x_init = -1, mean = 1)
# hist(test)

# test <- ars(f = dnorm, n = 1000, bounds = c(0,1), x_init = 0.5)
# hist(test, freq = F)
# curve(dnorm(x, 0, 1), 0, 1,  add = TRUE, col = "red")
# # 
# # 
test <- ars(f = dunif, n = 1000, bounds = c(10, 15), x_init = 11, min=10, max=15)
hist(test, freq = F)
curve(dunif(x, 10, 15), 10, 15,  add = TRUE, col = "red")

test <- ars(f = dexp, n = 1000, bounds = c(1, Inf), x_init = 1.5)
hist(test, freq = F)
curve(dexp(x, rate = 1/5), 0, 10,  add = TRUE, col = "red")s

#
test <- ars(f = dbeta, n = 1000, bounds = c(0, 5), x_init = 0.5, shape1 = 3, shape2 = 4)
hist(test, freq = F)
curve(dbeta(x, 3, 4), 0.01, .99, add = TRUE, col = "red")


hist(ars(dlaplace, 1000, x_init = -2, bounds = c(-5,-1)), s = 3)














test_that("check if the input density is function", {
  ## std. normal dist
  expect_equal(length(ars(dnorm, 100)), 100)

  ## exp. dist
  expect_equal(length(ars(dexp, 100, x_init = 0.5, c(0,Inf))), 100)

  ## gamma dist
  expect_equal(length(ars(dgamma, 100, x_init=10, bounds=c(0.0, Inf), shape=3.0, rate=2.0)), 100)

  ## case when the input density is not a function
  expect_error(ars(1, 100))
})

test_that("check if the bound is of length 2", {
  ## length 2 case
  expect_equal(length(ars(dnorm, 100, bounds = c(-100,100))), 100)

  ## case when length not equal to 2
  expect_error(ars(dnorm, 100, bounds = c(-1,0,1)))
})

test_that("check if the upper bound and the lower bound are equal", {
  #expect_warning(ars(dnorm,100,bounds = c(100,100)))
})

test_that("check if x_0 is between bounds", {
  ## when x_0 is at the left side of bounds
  expect_error(ars(dnorm, 100, x_init = -1, bounds = c(0,1)))

  ## when x_0 is at the right side of bounds
  expect_error(ars(dnorm, 100, x_init = 2, bounds = c(0,1)))

  ## when x_0 is in between bounds
  expect_equal(length(ars(dnorm, 100, x_init = 0.5, bounds = c(0,1))), 100)

})


test_that("check if the sampling distribution is close enough to the original distribution", {
  # Note we perform Kolmogorov-Smirnov Test with the null
  # hypothesis that the samples are drawn from the same
  # continuous distribution

  ## std. normal
  expect_equal(ks.test(ars(dnorm, 1000), rnorm(1000))$p.value > 0.05, T)

  ## exp(1)
  expect_equal(ks.test(ars(dexp, 1000, x_init = 5, bounds = c(0, Inf)), rexp(1000))$p.value > 0.05, T)

  ## gamma(3,2)
  expect_equal(ks.test(ars(dgamma, 1000, x_init = 5, bounds = c(0, Inf), shape = 3, scale = 2),
                       rgamma(1000, shape = 3, scale = 2))$p.value > 0.05, T)
  ## unif(0,1)
  expect_equal(ks.test(ars(dunif, 1000, x_init = 0.5, bounds = c(0,1)), runif(1000))$p.value > 0.05, T)

  ## logistics
  expect_equal(ks.test(ars(dlogis, 1000, x_init = 0, bounds = c(-10,10)), rlogis(1000))$p.value > 0.05, T)

  ## beta(3,2)
  expect_equal(ks.test(ars(dbeta, 100, x_init = 0.5, bounds = c(0, 1), shape1 = 3, shape2 = 2),
                       rbeta(100, shape1 = 3, shape2 = 2))$p.value > 0.05, T)

  # ## laplace
  # library(rmutil)
  # expect_equal(ks.test(ars(dlaplace, 1000, x_init = -2, bounds = c(-5,-1)),
  #                      rlaplace(1000))$p.value > 0.05, T)

  ## chi(2)
  expect_equal(ks.test(ars(dchisq, 100, x_init = 1, bounds = c(0, Inf), df = 2),
                       rchisq(100, df = 2))$p.value > 0.05, T)

  ## weibull(2,1)
  expect_equal(ks.test(ars(dweibull, 100, shape = 2, x_init = 1, bounds = c(0, Inf)),
                       rweibull(100, shape = 2))$p.value > 0.05, T)
})

test_that("check for non-log-concavity", {

  ## simple exponential exp(x^2)
  de = function (x) {
    return (exp(x^2))
  }
  expect_error(ars(de, 1000, x_init = 0, bounds = c(-5,5)))

  ## student t(2)
  expect_error(ars(dt, 1000, x_init = 1, bounds = c(-5,5), df = 2))

  ## cauchy
  expect_error(ars(dcauchy, 1000, x_init = 0, bounds = c(-5,5)))

  # pareto(1,2)
  expect_error(ars(dpareto, 1000, x_init = 3, bounds = c(1, Inf), m = 1, s = 2))

  ## lognormal
  expect_error(ars(dlnorm, 1000, x_init = 1, bounds = c(0, Inf)))
  #
  # ## F dist (1,1)
  expect_error(ars(stats::df, 1000, x_init = 1, bounds = c(0, Inf), df1 = 1, df2 = 2))
})

