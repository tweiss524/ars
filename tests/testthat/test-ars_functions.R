
# test for log_concave()
test_that('log_concave errors for non-decreasing vector', {

  # create a non-log-concave vector, i.e., non-decreasing
  non_concave1 <- c(3, 1, 5, -1, 0)
  expect_error(check_log_concave(non_concave1), "Function is not log-concave")

  non_concave2 <- c(0, 0, 0, 0.1)
  expect_error(check_log_concave(non_concave2), "Function is not log-concave")

  # constant vector should still be log-concave
  concave1 <- c(3, 3, 3, 3)
  expect_no_error(check_log_concave(concave1))

})



# test for initialize_abscissae()
test_that('initialize_abscissae returns sequence of proper length', {

  # case of infinite left bound
  # function will run loop until xinit > 0,
  # then output a vector of length 20
  xinit_1 <- -6
  bounds_1 <- c(-Inf, 0)
  h_1 <- function(x) {-x - 20}
  expect_equal(length(initialize_abscissae(xinit_1, h_1, bounds_1)), 20)


  # case of infinite right bound
  # function will run loop until xinit < 0,
  # then output a vector of length 20
  xinit_2 <- -6
  bounds_2 <- c(0, Inf)
  h_2 <- function(x) {-x}
  expect_equal(length(initialize_abscissae(xinit_2, h_2, bounds_2)), 20)

  # case of finite bounds
  xinit_3 <- -6
  bounds_3 <- c(-5, 5)
  h_3 <- function(x) {-2*x}
  expect_equal(length(initialize_abscissae(xinit_3, h_3, bounds_2)), 20)

})


# test for calc_z() function
test_that('test that calc_z computes calculation correctly', {

  # case when there is only 1 value in Tk, where
  # calc_z should return empty vector
  tk_1 <- 1
  h_tk_1 <- -3.751
  hprime_tk_1 <- 0.069

  calc_z_with_1_value <- calc_z(tk_1, h_Tk_1, hprime_Tk_1)
  expect_equal(calc_z_with_1_value, c())

  # case where Tk has more than 1 value
  tk_2 <- c(0.000, 1.052)
  h_tk_2 <- c(3.751, -3.682)
  hprime_tk_2 <- c(0.069, 0.062)

  # manual calculation
  z_man <- (-3.682 - 3.751 - 1.052*0.062 + 0.000*0.069) / (0.069-0.062)

  # function calculation
  z_func <- calc_z(tk_2, h_tk_2, hprime_tk_2)

  expect_equal(z_man, z_func)

})


# test for u() function
test_that('test that u computes calculations correctly', {

  tk <- c(0.000, 1.052)
  h_tk <- c(-3.751, -3.682)
  hprime_tk <- c(0.069, 0.062)
  x_1 <- 0.3

  # case when zk is an empty vector
  zk_1 <- c()

  # manual calculation
  u_man_1 <- -3.751 + (x_1 - 0.000) * 0.069

  # function calculation
  u_func_1 <- u(x_1, zk_1, tk, h_tk, hprime_tk)

  expect_equal(u_man_1, u_func_1)


  # 2 tests for case when zk is not empty:
  zk_2 <- 0.526

  # case when x is on the left hand side of zk
  x_left <- 0.1
  u_man_left <- -3.751 + (x_left - 0.000)*0.069
  u_func_left <- u(x_left, zk_2, tk, h_tk, hprime_tk)
  expect_equal(u_man_left, u_func_left)

  # case when x is on the right hand side of zk
  x_right <- 0.6
  u_man_right <- -3.682 + (x_right - 1.052) * 0.062
  u_func_right <- u(x_right, zk_2, tk, h_tk, hprime_tk)
  expect_equal(u_man_right, u_func_right)

})


# test for l() function
test_that('test that l computes calculations correctly',{

  # case when Tk only has one element.

  x_1 <- 0.3
  zk_1 <- NULL
  tk_1 <- 0.000
  h_tk_1 <- -3.751
  hprime_tk_1 <- 0.069

  # manual calculation
  l_man_1 <- -3.751 + (x_1 - 0.000)*0.069

  # function calculation
  l_func_1 <- l(x_1, zk_1, tk_1, h_tk_1, hprime_tk_1)

  expect_equal(l_man_1, l_func_1)


  # case when x < Tk[1]
  # should return -Inf
  x_less <- -1
  zk_less <- 0.562
  tk_less <- c(0.000, 1.052)
  h_tk_less <- c(-3.751, -3.682)
  hprime_tk_less <- c(0.069, 0.062)

  l_less <- -Inf
  l_func_less <- l(x_less, zk_less, tk_less, h_tk_less, hprime_tk_less)
  expect_equal(l_less, l_func_less)

  # case when x > Tk[k]
  # should return -Inf
  x_gre <- 1.6
  zk_gre <- 0.414
  tk_gre <- c(0.012, 1.09)
  h_tk_gre <- c(-3.99, -3.518)
  hprime_tk_gre <- c(0.071, 0.094)

  l_gre <- -Inf
  l_func_gre <- l(x_gre, zk_gre, tk_gre, h_tk_gre, hprime_tk_gre)
  expect_equal(l_gre,l_func_gre)


  # case when x is on the left hand side of zk but bigger than Tk[1]

  x_l <- 0.4
  zk_l <- 0.562
  tk_l <- c(0.000, 1.052)
  h_tk_l <- c(-3.751, -3.682)
  hprime_tk_l <- c(0.069, 0.062)

  l_man_l <- ((1.052 - x_l)*(-3.751) + (x_l - 0.000)*(-3.682)) / (1.052 - 0.000)
  l_func_l <- l(x_l, zk_l, tk_l, h_tk_l, hprime_tk_l)
  expect_equal(l_man_l, l_func_l)


  # case when x in [Tk_j, Tk_j+1]

  x_in <- 1.82
  zk_in <- c(0.562, 0.333)
  tk_in <- c(0.000, 1.052, 2.105)
  h_tk_in <- c(-3.751, -3.682, -3.518)
  hprime_tk_in <- c(0.069, 0.062, 0.053)

  l_man_in <- ((2.105 - x_in)*(-3.682) + (x_in - 1.052)*(-3.518)) / (2.105-1.052)
  l_func_in <- l(x_in, zk_in, tk_in, h_tk_in, hprime_tk_in)
  expect_equal(l_man_in, l_func_in)

})


# test for calc_probs()
test_that('test the calc_probs function computes calcuations correctly',{


  tk_1 <- c(0.000,1.052)
  zk_1 <- 0.562
  h_tk_1 <- c(-3.751,-3.682)

  ## The first element
  hprime_tk_1 <- c(0, 0.062)
  bounds_1 <- c(0, 4)

  # simulating calculations
  z_all_1 <- c(bounds_1[1], zk_1, bounds_1[2])
  num_bins_1 <- length(z_all_1) - 1
  z_1_1 <- z_all_1[1:num_bins_1]
  z_2_1 <- z_all_1[2:length(z_all_1)]
  u_z1_1 <- u(z_1_1, zk_1, tk_1, h_tk_1, hprime_tk_1) # c(-3.751,-3.71238)
  u_z2_1 <- u(z_2_1, zk_1, tk_1, h_tk_1, hprime_tk_1) # c(-3.71238,-3.499224)

  # manual
  unnormalized_prob_man_1 <- c(exp(-3.751)*(0.562-0),
                               (exp(-3.499224) - exp(-3.71238))/0.062)

  # function
  unnormalized_prob_fun_1 <- calc_probs(tk_1, zk_1, h_tk_1, hprime_tk_1, bounds_1)[[2]]


  normalized_prob_man_1 <- unnormalized_prob_man_1/sum(unnormalized_prob_man_1)
  normalized_prob_fun_1 <- calc_probs(tk_1, zk_1, h_tk_1, hprime_tk_1, bounds_1)[[1]]


  expect_equal(normalized_prob_man_1[1],normalized_prob_fun_1[1])
  expect_equal(normalized_prob_man_1[2],normalized_prob_fun_1[2])
  expect_equal(unnormalized_prob_man_1[1],unnormalized_prob_fun_1[1])
  expect_equal(unnormalized_prob_man_1[2],unnormalized_prob_fun_1[2])


})



# test for sample_sk() function
test_that('test that sample_sk computes calculations correctly',{


  tk_1 <- c(0.000, 1.052)
  zk_1 <- c(0.562)
  h_tk_1 <- c(-3.751, -3.682)

  ## The first element
  hprime_tk_1 <- c(0.062, 0)
  bounds_1 <- c(0, 4)

  z_all_1 <- c(bounds_1[1], zk_1, bounds_1[2])
  num_bins_1 <- length(z_all_1) - 1
  z_1_1 <- z_all_1[1:num_bins_1]
  u_z1_1 <- u(z_1_1, zk_1, tk_1, h_tk_1, hprime_tk_1) # c(-3.751,-3.682)
  p <- calc_probs(tk_1, zk_1, h_tk_1, hprime_tk_1, bounds_1)

  # normalized
  prob <- p[[1]] ## c(0.2382642, 0.7617358)

  # unnormalized
  unnorm_prob <- p[[2]] ## c(0.02706999, 0.08654333)


  unif <- 0.937642
  i <- 2

  # manual computation
  xstar_man_2 <- 0.562 + 0.937642*(4-0.562)
  set.seed(6)

  # function
  xstar_fun_2 <- sample_sk(tk_1, zk_1, h_tk_1, hprime_tk_1, bounds_1)

  # test that xstars are the same
  expect_equal(sum(which(abs(xstar_man_2 - xstar_fun_2) < 1e-7)), 1)

})
