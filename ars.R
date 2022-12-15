
#' Adaptive Rejection Sampling
#'
#' \code{ars}
#'
#' @description This function implements an adaptive rejection sampling algorithm to generate n samples from any univariate log-concave function.
#'
#' @author Jeffrey Kuo, Hsiang-Chuan Sha, Tessa Weiss
#'
#' @usage ars(f, n, bounds = c(-Inf, Inf), x_init = NA)
#'
#'
#' @param f function from which to generate samples from
#' @param n numbers of desired samples
#' @param bounds vector of length 2 defining the lower and upper bound of of the distribution f
#' @param x_init initial point within the bounds that used to initialize the abscissae points
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
#' ars(dnorm, n = 10, bounds = c(-10,10), x_init = 1)
#'
#' [1] -0.2133253 -2.4711098 -0.3627347 -0.2517910 -1.8091649  0.8678341 -2.0045281 -1.1447111 -2.1649985 -0.4544712
#'
#' ## plot 5000 points from a Gamma(3, 4) distribution
#' gam_samps <- ars(f = function(x) {dgamma(x, 3, 4)}, n = 5000, bounds = c(0, Inf), x_init = 1)
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

  source("ars_functions.R")

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

  Tk <- initialize_abcissae(x_init, hprime, bounds)
  # print("Tk:")
  # print(Tk)

  num_samps <- 1
  #print("Found absissae")
  h_Tk <- sapply(Tk, h)
  #print("Found h")
  # print(h_Tk)

  # check whether the function is defined inside bounds by removing infinite log(f(x))
  Tk <- Tk[is.finite(h_Tk)]
  h_Tk <- h_Tk[is.finite(h_Tk)]
  assertthat::assert_that(length(h_Tk) > 0, msg = "Function not defined in bounds")
  # print("newTk")
  # print(Tk)

  hprime_Tk <- sapply(Tk, hprime)
  # print("Found derivative")
  # print(hprime_Tk)

  # remove x such that h'(x) is infinite
  h_Tk <- h_Tk[!is.na(hprime_Tk)]
  Tk <- Tk[!is.na(hprime_Tk)]
  hprime_Tk <- hprime_Tk[!is.na(hprime_Tk)]

  # print("hprime_Tk:")
  # print(hprime_Tk)
  #print(sum(abs(hprime_Tk[2:length(hprime_Tk)] - hprime_Tk[1:(length(hprime_Tk) -1)]) <=  1e-8)  == (length(hprime_Tk)-1))


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
  # print("zk:")
  # print(zk)


  ########## SAMPLING STEP ##########

  while(num_samps <= n) {

    # sample xstar from sk
    xstar <- sample_sk(Tk, zk, h_Tk, hprime_Tk, bounds)
    #print(paste("xstar:", xstar))
    #print(paste("l:", l(xstar, Tk, h_Tk, hprime_Tk) ))
    #print(paste("h:", h(xstar)))
    # print(zk)
    # print(Tk)
    #print(paste("u", u(xstar, zk, Tk, h_Tk, hprime_Tk)))
    # print("zk")
    # print(zk)
    assertthat::assert_that((xstar >= bounds[1]) && (xstar <= bounds[2]), msg = "sampled xstar not in bounds")
    #zk <- sort(zk)
    #assertthat::assert_that((l(xstar, Tk, h_Tk, hprime_Tk) <= h(xstar)) && (h(xstar) <= u(xstar, zk, Tk, h_Tk, hprime_Tk)), msg = "lhu test: Not log concave")
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
      if(sum(abs(hprime_Tk[2:len_hptk_new] - hprime_Tk[1:(len_hptk_new-1)]) <=  1e-8)  == (len_hptk_new-2)){
        hprime_Tk <- as.integer(round(hprime_Tk))
        ind <- c(which(hprime_Tk == unique(hprime_Tk)[1])[1],which(hprime_Tk == unique(hprime_Tk)[2])[1])
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
test <- ars(f = function(x){dnorm(x, 0, 1)}, n = 1000, x_init = 1, bounds = c(-Inf,Inf))
# test
hist(test)
hist(ars(f = function(x){dexp(x, 4)}, n = 100, x_init = 1, bounds= c(0, Inf)), freq = F)
#
# getwd()
#
# # curve(dgamma(x, 3, 4), 0, add = TRUE, col = "red")
# # ars(f = dunif, n = 1000, bounds = c(0,1))
# #
# #
# #test <- ars(f = dnorm, n = 1000, mean = 10000, x_init = 10000)
#
# # hist(ars(dcauchy, 1000, x_init = 0, bounds = c(-5,5)))
# # hist(test, freq = F)
# # curve(dnorm(x, 50, 1), 40, 60,  add = TRUE, col = "red")
# #
# # #
# # # # testing with gamma
# # test <- ars(f = dgamma, n = 1000,  bounds = c(1,10), x_init = 5, shape = 9, rate = 2)
# # #
# # hist(test, freq = F)
# # curve(dgamma(x, 9, 2), 1, 10, add = TRUE, col = "red")
# #
# # #
# # # # testing with logistics
# #test <- ars(f = dlogis, n = 1000,  bounds = c(1,5), x_init = 2)
# #hist(test, freq = F)
# #curve(dlogis(x, 0, 1), 1, 5, add = TRUE, col = "red")
# #
# # #
# # test <- ars(f = dnorm, n = 1000,  bounds = c(-Inf,0), x_init = -1, mean = 1)
# # hist(test)
#
# # test <- ars(f = dnorm, n = 1000, bounds = c(0,1), x_init = 0.5)
# # hist(test, freq = F)
# # curve(dnorm(x, 0, 1), 0, 1,  add = TRUE, col = "red")
# # #
# # #
# test <- ars(f = function(x){dnorm(x, 5, 1)}, n = 1000, bounds = c(10, 15), x_init =11)
# hist(test, freq = F)
# curve(dunif(x, 10, 15), 10, 15,  add = TRUE, col = "red")
#
# test <- ars(f = dexp, n = 1000, bounds = c(1, Inf), x_init = 1.5)
# hist(test, freq = F)
# curve(dexp(x, rate = 1/5), 0, 10,  add = TRUE, col = "red")
#
# #
# test <- ars(f = dbeta, n = 1000, x_init = 0.5, shape1 = 4, shape2 = 3)
# hist(test, freq = F)
# curve(dbeta(x, 4, 3), 0.01, .99, add = TRUE, col = "red")
#
#
# hist(ars(dlaplace, 1000, x_init = -2))
# hist(ars(dlaplace, 1000, x_init = -2, bounds = c(-5,1)), s = 3)
#
#
# nor_test <- ars(dnorm, n = 1000)
# unif_test <- ars(dunif, n = 1000, bounds = c(10,15), x_init = 11, min = 10, max = 15)
# laplace_test <- ars(dlaplace, n = 1000, x_init = -2)
# beta_test <- ars(dbeta, bounds = c(0, 1), x_init = 0.5, shape1 = 4, shape2 = 3)
# exp_test <- ars(dexp, bounds = c(0, Inf), x_init = 4, rate = 5)
#
#
# abs(mean(rnorm(1000)) - mean(nor_test)) < 0.01
# abs(var(rnorm(1000)) - var(nor_test)) < 0.01
#
#
#
#
#
#
#
#
#
#
#
# test_that("provides samples close to actual distribution", {
#
#   # implements KS test
#
#   # N(0, 1)
#   nor_ars <- ars(dnorm, n = 1000)
#   true_nor <- rnorm(1000, 0, 1)
#   expect_equal(ks.test(nor_ars, true_nor)$p.value <= 0.05, FALSE)
#
#
#   # Unif(10, 15)
#   f_unif <- function(x) {return(dunif(x, 10, 15))}
#   unif_ars <- ars(f_unif, n = 1000, bounds = c(10,15), x_init = 11)
#   true_unif <- runif(1000, 10, 15)
#   expect_equal(ks.test(unif_ars, true_unif)$p.value <= 0.05, FALSE)
#
#   # Beta(4, 3)
#   f_beta <- function(x) {return(dbeta(x, 4, 3))}
#   beta_ars <- ars(f_beta, n = 1000, bounds = c(0, 1), x_init = 0.5)
#   true_beta <- rbeta(1000, 4, 3)
#   expect_equal(ks.test(beta_ars, true_beta)$p.value <= 0.05, FALSE)
#
#   # Exp(5)
#   f_exp <- function(x) {return(dexp(x, 5))}
#   exp_ars <- ars(f_exp, n = 1000, bounds = c(0, Inf), x_init = 4)
#   true_exp <- rexp(1000, 5)
#   expect_equal(ks.test(exp_ars, true_exp)$p.value <= 0.05, FALSE)
#
#   # Chi-sqared(4)
#   f_chi <- function(x) {return(dchisq(x, 4))}
#   chi_ars <- ars(f_chi, n = 1000, bounds = c(0, Inf), x_init = 1)
#   true_chi <- rchisq(1000, 4)
#   expect_equal(ks.test(chi_ars, true_chi)$p.value <= 0.05, FALSE)
#
# })
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
# test_that("check if the input density is function", {
#   ## std. normal dist
#   expect_equal(length(ars(dnorm, 100)), 100)
#
#   ## exp. dist
#   expect_equal(length(ars(dexp, 100, x_init = 0.5, c(0,Inf))), 100)
#
#   ## gamma dist
#   expect_equal(length(ars(dgamma, 100, x_init=10, bounds=c(0.0, Inf), shape=3.0, rate=2.0)), 100)
#
#   ## case when the input density is not a function
#   expect_error(ars(1, 100))
# })
#
# test_that("check if the bound is of length 2", {
#   ## length 2 case
#   expect_equal(length(ars(dnorm, 100, bounds = c(-100,100))), 100)
#
#   ## case when length not equal to 2
#   expect_error(ars(dnorm, 100, bounds = c(-1,0,1)))
# })
#
# test_that("check if the upper bound and the lower bound are equal", {
#   expect_warning(ars(dnorm,100,bounds = c(100,100)))
# })
#
# test_that("check if x_0 is between bounds", {
#   ## when x_0 is at the left side of bounds
#   expect_error(ars(dnorm, 100, x_init = -1, bounds = c(0,1)))
#
#   ## when x_0 is at the right side of bounds
#   expect_error(ars(dnorm, 100, x_init = 2, bounds = c(0,1)))
#
#   ## when x_0 is in between bounds
#   expect_equal(length(ars(dnorm, 100, x_init = 0.5, bounds = c(0,1))), 100)
#
# })
#
#
# test_that("check if the sampling distribution is close enough to the original distribution", {
#   # Note we perform Kolmogorov-Smirnov Test with the null
#   # hypothesis that the samples are drawn from the same
#   # continuous distribution
#
#   ## std. normal
#   expect_equal(ks.test(ars(dnorm, 1000), rnorm(1000))$p.value > 0.05, T)
#
#   ## exp(1)
#   expect_equal(ks.test(ars(dexp, 1000, x_init = 5, bounds = c(0, Inf)), rexp(1000))$p.value > 0.05, T)
#
#   ## gamma(3,2)
#   f <- function(x){return(dgamma(x, 3, 2))}
#
#   expect_equal(ks.test(ars(f, 1000, x_init = 5, bounds = c(0, Inf)),
#                        rgamma(1000, shape = 3, scale = 2))$p.value > 0.05, T)
#   ## unif(0,1)
#   expect_equal(ks.test(ars(dunif, 1000, x_init = 0.5, bounds = c(0,1)), runif(1000))$p.value > 0.05, T)
#
#   ## logistics
#   expect_equal(ks.test(ars(dlogis, 1000, x_init = 0, bounds = c(-10,10)), rlogis(1000))$p.value > 0.05, T)
#
#   ## beta(3,2)
#   f <- function(x){return(dbeta(x, 3, 2))}
#   expect_equal(ks.test(ars(f, 100, bounds = c(0, 1)),
#                        rbeta(100, shape1 = 3, shape2 = 2))$p.value > 0.05, T)
#
#   # ## laplace
#   library(rmutil)
#   expect_equal(ks.test(ars(dlaplace, 1000, x_init = -2, bounds = c(-5,5)),
#                        rlaplace(1000))$p.value > 0.05, T)
#
#   ## chi(2)
#   f <- function(x){return(dchisq(x, 2))}
#   expect_equal(ks.test(ars(f, 100, x_init = 1, bounds = c(0, Inf)),
#                        rchisq(100, df = 2))$p.value > 0.05, T)
#
#   ## weibull(2,1)
#   f <- function(x){return(dweibull(x, 2))}
#   expect_equal(ks.test(ars(f, 100, x_init = 1, bounds = c(0, Inf)),
#                        rweibull(100, shape = 2))$p.value > 0.05, T)
# })
#
# test_that("check for non-log-concavity", {
#
#   ## simple exponential exp(x^2)
#   de = function (x) {
#     return (exp(x^2))
#   }
#   expect_error(ars(de, 1000, x_init = 0, bounds = c(-5,5)))
#
#   ## student t(2)
#   expect_error(ars(dt, 1000, x_init = 1, bounds = c(-5,5), df = 2))
#
#   ## cauchy
#   expect_error(ars(dcauchy, 1000, x_init = 0, bounds = c(-5,5)))
#
#   # pareto(1,2)
#   expect_error(ars(dpareto, 1000, x_init = 3, bounds = c(1, Inf), m = 1, s = 2))
#
#   ## lognormal
#   expect_error(ars(dlnorm, 1000, x_init = 1, bounds = c(0, Inf)))
#   #
#   # ## F dist (1,1)
#   expect_error(ars(stats::df, 1000, x_init = 1, bounds = c(0, Inf), df1 = 1, df2 = 2))
# })
#
#
#
#
# test_that("function errors out for non log-concave functions", {
#
#   # cauchy dist
#   expect_error(ars(dcauchy, n = 1000, x_init = 1, bounds = c(-10,10)), "Function is not log-concave")
#
#   # t-dist
#   expect_error(ars(dt, 1000, bounds = c(-5,5), x_init = 1, df = 4),"Function is not log-concave")
#
#   # pareto dist
#   expect_error(ars(dpareto, n = 1000, bounds = c(1, Inf), x_init = 2, m = 6, s = 2), "Function is not log-concave")
#
#   # lognormal dist
#   expect_error(ars(dlnorm, n = 1000, bounds = c(0, Inf), x_init = 1), "Function is not log-concave")
# })
#
# test_that("function detects invalid inputs", {
#
#   # f is not a function
#   expect_error(ars(f = "hi", n = 1000, x_init = 1, bounds = c(-10,10)), "f must be a function")
#
#   # n is not numeric
#   expect_error(ars(dnorm, n = "1000", bounds = c(-5,5), x_init = 1), "n must be numeric")
#
#   # bounds not numeric
#   expect_error(ars(dnorm, n = 1000, bounds = c("hi", "hello"), x_init = 2), "Bounds must be numeric vector of length 2")
#
#   # bounds are not in order
#   expect_warning(ars(dnorm, n = 1000, bounds = c(Inf, 1), x_init = 2))
#
#   # x_init not in bounds
#   expect_error(ars(dnorm, n = 1000, bounds = c(0, Inf), x_init = -5), "x_init point must be inside bounds")
# })
#
#
# f <- function(x){dbeta(x, 10, 5)}
# optim(1, f, method = 'L-BFGS-B', control = list(fnscale=-1))


