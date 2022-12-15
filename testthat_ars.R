


library(testthat)


# Functions for testing are helper functions within ars:

# initialize abscsissae
# u, l, z, calc_probs, sample_sk



test_that("errors out for non log-concave functions", {
  
  # cauchy dist
  f_cauchy <- function(x){dcauchy(x, 2,4)}
  expect_error(ars(f_cauchy, n = 1000, x_init = 1, bounds = c(-10,10)), "Function is not log-concave")
  
  # t-dist
  f_t <- function(x) {dt(x, 4)}
  expect_error(ars(f_t, 1000, bounds = c(-5,5), x_init = 1), "Function is not log-concave")
  
  # pareto dist
  f_par <- function(x) {rmutil::dpareto(x, 6, 2)}
  expect_error(ars(f_par, n = 1000, bounds = c(1, Inf), x_init = 2), "Function is not log-concave")
  
  # lognormal dist
  expect_error(ars(dlnorm, n = 1000, bounds = c(0, Inf), x_init = 1), "Function is not log-concave")
})

test_that("detects invalid inputs", {
  
  # f is not a function 
  expect_error(ars(f = "hi", n = 1000, x_init = 1, bounds = c(-10,10)), "f must be a function")
  
  # n is not numeric
  expect_error(ars(dnorm, n = "1000", bounds = c(-5,5), x_init = 1), "n must be numeric")
  
  # bounds not numeric
  expect_error(ars(dnorm, n = 1000, bounds = c("hi", "hello"), x_init = 2), "Bounds must be numeric vector of length 2")
  
  # bounds are not in order
  expect_warning(ars(dnorm, n = 1000, bounds = c(Inf, 1), x_init = 2))
  
  # x_init not in bounds
  expect_error(ars(dnorm, n = 1000, bounds = c(0, Inf), x_init = -5), "x_init point must be inside bounds")
})


test_that("provides samples close to actual distribution", {
  
  # implements KS test
  
  # N(0, 1)
  nor_ars <- ars(dnorm, n = 1000)
  true_nor <- rnorm(1000, 0, 1)
  expect_equal(ks.test(nor_ars, true_nor)$p.value <= 0.05, FALSE)
  
  
  # Unif(10, 15)
  f_unif <- function(x) {return(dunif(x, 10, 15))}
  unif_ars <- ars(f_unif, n = 1000, bounds = c(10,15), x_init = 11)
  true_unif <- runif(1000, 10, 15)
  expect_equal(ks.test(unif_ars, true_unif)$p.value <= 0.05, FALSE)
  
  # Beta(4, 3)
  f_beta <- function(x) {return(dbeta(x, 4, 3))}
  beta_ars <- ars(f_beta, n = 1000, bounds = c(0, 1), x_init = 0.5)
  true_beta <- rbeta(1000, 4, 3)
  expect_equal(ks.test(beta_ars, true_beta)$p.value <= 0.05, FALSE)
  
  # Exp(5)
  f_exp <- function(x) {return(dexp(x, 5))}
  exp_ars <- ars(f_exp, n = 1000, bounds = c(0, Inf), x_init = 4)
  true_exp <- rexp(1000, 5)
  expect_equal(ks.test(exp_ars, true_exp)$p.value <= 0.05, FALSE)
  
  # Chi-sqared(4)
  f_chi <- function(x) {return(dchisq(x, 4))}
  chi_ars <- ars(f_chi, n = 1000, bounds = c(0, Inf), x_init = 1)
  true_chi <- rchisq(1000, 4)
  expect_equal(ks.test(chi_ars, true_chi)$p.value <= 0.05, FALSE)
  
})
