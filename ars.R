
##### MAIN ARS FUNCTION #####
#############################



#' Adaptive Rejection Sampling
#'
#' \code{ars}
#'
#' @description This function implements an adaptive rejection sampling algorithm to generate a vector of samples from any univariate log-concave function.
#'
#' @author Jeffrey Kuo, Hsiang-Chuan Sha, Tessa Weiss
#'
#' @usage ars(f, n = 1000, bounds = c(-Inf, Inf), x_init = NA)
#'
#'
#' @param f function from which to generate samples from
#' @param n numbers of desired samples; default is 1000
#' @param bounds vector of length 2 defining the lower and upper bounds of the distribution f; default is (-Inf, Inf)
#' @param x_init initial point within the bounds that is used to initialize the abscissae points; default is NA
#'
#' @return returns a vector of n samples from the distribution function f
#'
#'
#' @references
#' Gilks, W. R., & Wild, P. (1992). Adaptive rejection sampling for Gibbs sampling.
#' \emph{Journal of the Royal Statistical Society: Series C (Applied Statistics)}, 41(2), 337-348.
#'
#' @examples
#' ## sample 1000 points from a N(0, 1) distribution
#' ars(dnorm, n = 1000, bounds = c(-10, 10), x_init = 1)
#'
#'
#' ## plot 5000 points from a Gamma(3, 4) distribution
#' func <- function(x) {dgamma(x, 3, 4)}
#' gam_samps <- ars(f = func, n = 5000, bounds = c(0, Inf), x_init = 1)
#' hist(gam_samps)
#'
#' @importFrom stats runif
#' @importFrom stats optim
#' @importFrom numDeriv grad
#' @importFrom assertthat assert_that
#' @importFrom rmutil dpareto
#'
#' @export
#'



ars <- function(f, n = 1000, bounds = c(-Inf, Inf), x_init = NA) {

  #source("ars_functions.R")

  # if no initial point provided, take the mode of the function f as x_init
  if (is.na(x_init)) {
    x_init <- optim(1, f, method = 'L-BFGS-B', control = list(fnscale=-1))$par
  }

  # assertions on input parameters
  assertthat::assert_that(is.function(f), msg = "f must be a function")
  assertthat::assert_that((is.numeric(n) & (n > 0)), msg = "n must be a positive integer")
  assertthat::assert_that(is.numeric(x_init), msg = "x_init must be numeric")
  assertthat::assert_that(f(x_init) > 1e-8, msg = "x_init must have positive probability at f")
  assertthat::assert_that(is.vector(bounds) && (length(bounds) == 2) && (is.numeric(bounds)), msg = "Bounds must be numeric vector of length 2")

  # automatically reverse the bounds if lower bound is greater than upper bound

  if (bounds[1] > bounds[2]) {

    bounds <- c(bounds[2], bounds[1])

    warning(paste0("Lower bound must be smaller than upper bound. ",
                   sprintf("New bounds are (%s, %s)", bounds[1], bounds[2])))

  }


  # automatically set the bounds to the default (-Inf, Inf) if they are the same value

  if (bounds[1] == bounds[2]) {

    bounds <- c(-Inf, Inf)

    warning(paste0("Lower bound must be smaller than upper bound. ",
                   sprintf("New bounds are (%s, %s)", bounds[1], bounds[2])))

  }

  assertthat::assert_that((x_init >= bounds[1]) && (x_init <= bounds[2]), msg = "x_init point must be inside bounds")

  # make sure n is an integer
  n <- as.integer(n)


  # function h takes x as input and outputs log(f(x))
  h <- function(x) {

    logf <- log(f(x))
    #assertthat::assert_that(is.finite(logf), msg = sprintf("Invalid Bounds case logf: %s.", x))
    return(logf)
  }

  # function hprime takes x as input and outputs h'(x)
  hprime <- function(x) {
    der <- numDeriv::grad(h, x)
    #der <- (h(x + 1e-8) - h(x))/1e-8
    return(der)
  }

  # initializing sample vector
  samps <- rep(NA, n)


  ########## INITIALIZING STEP ##########

  Tk <- initialize_abscissae(x_init, hprime, bounds)

  num_samps <- 1
  
  h_Tk <- sapply(Tk, h)

  # check whether the function is defined inside bounds by removing infinite log(f(x))
  Tk <- Tk[is.finite(h_Tk)]
  h_Tk <- h_Tk[is.finite(h_Tk)]
  assertthat::assert_that(length(h_Tk) > 0, msg = "Function not defined in bounds")

  # derivatives at abscissae
  hprime_Tk <- sapply(Tk, hprime)

  # remove x such that h'(x) is infinite
  h_Tk <- h_Tk[!is.na(hprime_Tk)]
  Tk <- Tk[!is.na(hprime_Tk)]
  hprime_Tk <- hprime_Tk[!is.na(hprime_Tk)]

  # handle the case when hprime_Tk are all equal to a constant (same value)
  len_hptk <- length(hprime_Tk)

  if (sum(abs(hprime_Tk[2:len_hptk] - hprime_Tk[1:(len_hptk-1)]) <=  1e-8)  == (len_hptk - 1)) {

    # keep the first non-NA element
    hprime_Tk <- hprime_Tk[2]
    Tk <- Tk[2]
    h_Tk <- h_Tk[2]

  } else if (sum(abs(hprime_Tk[2:len_hptk] - hprime_Tk[1:(len_hptk-1)]) <=  1e-8)  == (len_hptk - 2)){

    # keep the first element for each repetitive hprime_Tk if hprime_Tk only takes values of two constants
    hprime_Tk <- as.integer(round(hprime_Tk))
    ind <- c(which(hprime_Tk == unique(hprime_Tk)[1])[1], which(hprime_Tk == unique(hprime_Tk)[2])[1])
    hprime_Tk <- hprime_Tk[ind]
    Tk <- Tk[ind]
    h_Tk <- h_Tk[ind]
  }

  check_log_concave(hprime_Tk)


  zk <- calc_z(Tk, h_Tk, hprime_Tk)


  ########## SAMPLING STEP ##########

  while(num_samps <= n) {

    # sample xstar from sk
    xstar <- sample_sk(Tk, zk, h_Tk, hprime_Tk, bounds)
    assertthat::assert_that((xstar >= bounds[1]) && (xstar <= bounds[2]), msg = "sampled xstar not in bounds")
    w <- runif(1)

    # squeezing test
    if(w <= exp(l(xstar, zk, Tk, h_Tk, hprime_Tk) - u(xstar, zk, Tk, h_Tk, hprime_Tk))) {

      samps[num_samps] <- xstar
      num_samps <- num_samps + 1
    }

    # rejection test
    else {
      
      h_xstar <- h(xstar)
      if(!is.finite(h_xstar)) { break }
      hprime_xstar <- hprime(xstar)


      ########## UPDATING STEP ##########

      if (w <= exp(h_xstar - u(xstar, zk, Tk, h_Tk, hprime_Tk))) {
        
        samps[num_samps] <- xstar
        num_samps <- num_samps + 1
      }

      # append the x_star that fails the squeezing test into Tk,
      # and calculate the corresponding h_Tk+1, hprime_Tk+1
      Tk <- sort(c(Tk, xstar))
      h_Tk <- append(h_Tk, h_xstar, after = (which(Tk == xstar) - 1))
      hprime_Tk <- append(hprime_Tk, hprime_xstar, after = (which(Tk == xstar) - 1))

      h_Tk <- h_Tk[!is.na(hprime_Tk)]
      Tk <- Tk[!is.na(hprime_Tk)]
      hprime_Tk <- hprime_Tk[!is.na(hprime_Tk)]

      # check the behavior of hprime_Tk
      len_hptk_new <- length(hprime_Tk)
      if(sum(abs(hprime_Tk[2:len_hptk_new] - hprime_Tk[1:(len_hptk_new-1)]) <=  1e-8)  == (len_hptk_new-2)) {
        
        hprime_Tk <- as.integer(round(hprime_Tk))
        ind <- c(which(hprime_Tk == unique(hprime_Tk)[1])[1], which(hprime_Tk == unique(hprime_Tk)[2])[1])
        hprime_Tk <- hprime_Tk[ind]
        Tk <- Tk[ind]
        h_Tk <- h_Tk[ind]
        
      }

      check_log_concave(hprime_Tk)

      # update the zk+1
      zk <- calc_z(Tk, h_Tk, hprime_Tk)

    } # end else
  } # end while

  return(samps)
}


#
#
# # # testing with normal
# #
#test <- ars(f = function(x){dnorm(x, 0, 1)}, n = 1000, x_init = 1, bounds = c(-Inf,Inf))
# test
#hist(test)
#hist(ars(f = function(x){dexp(x, 4)}, n = 100, x_init = 1, bounds= c(0, Inf)), freq = F)
#
# getwd()
#
# # curve(dgamma(x, 3, 4), 0, add = TRUE, col = "red")


