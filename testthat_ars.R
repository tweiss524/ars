


library(test_that)


# Functions for testing are helper functions within ars:

# initialize abscsissae
# u, l, z, calc_probs, sample_sk



test_that("errors out for non log-concave functions", {
  
  # cauchy dist
  expect_error(ars(dcauchy, n = 1000, x_init = 1, bounds = c(-10,10)), "Function is not log-concave")
  
  # t-dist
  expect_error(ars(dt, 1000, bounds = c(-5,5), x_init = 1, df = 4), "Function is not log-concave")
  
  # pareto dist
  expect_error(ars(dpareto, n = 1000, bounds = c(1, Inf), x_init = 2, m = 6, s = 2), "Function is not log-concave")
  
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
  unif_ars <- ars(dunif, n = 1000, bounds = c(10,15), x_init = 11, min = 10, max = 15)
  true_unif <- runif(1000, 10, 15)
  expect_equal(ks.test(unif_ars, true_unif)$p.value <= 0.05, FALSE)
  
  # Beta(4, 3)
  beta_ars <- ars(dbeta, n = 1000, bounds = c(0, 1), x_init = 0.5, shape1 = 4, shape2 = 3)
  true_beta <- rbeta(1000, 4, 3)
  expect_equal(ks.test(beta_ars, true_beta)$p.value <= 0.05, FALSE)
  
  # Exp(5)
  exp_ars <- ars(dexp, n = 1000, bounds = c(0, Inf), x_init = 4, rate = 5)
  true_exp <- rexp(1000, 5)
  expect_equal(ks.test(exp_ars, true_exp)$p.value <= 0.05, FALSE)
  
  # Chi-squared(4)
  chi_ars <- ars(dchisq, n = 1000, bounds = c(0, Inf), x_init = 1, df = 4)
  true_chi <- rchisq(1000, 4)
  expect_equal(ks.test(chi_ars, true_chi)$p.value <= 0.05, FALSE)
  
})
