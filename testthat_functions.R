
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
  expect_equal(length(initialize_abcissae(xinit_1, h_1, bounds_1)), 20)
  
  
  # case of infinite right bound
  # function will run loop until xinit < 0,
  # then output a vector of length 20
  xinit_2 <- -6
  bounds_2 <- c(0, Inf)
  h_2 <- function(x) {-x}
  expect_equal(length(initialize_abcissae(xinit_2, h_2, bounds_2)), 20)
  
  # case of finite bounds
  xinit_3 <- -6
  bounds_3 <- c(-5, 5)
  h_3 <- function(x) {-2*x}
  expect_equal(length(initialize_abcissae(xinit_3, h_3, bounds_2)), 20)
  
})


# test for calc_z()
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
  
  calc_z_normal_case <- calc_z(tk_2, h_tk_2, hprime_tk_2)
  
  # manual calculation
  z_val <- (-3.682 - 3.751 - 1.052*0.062 + 0.000*0.069) / (0.069-0.062)
  
  expect_equal(calc_z_normal_case, z_val)
  
})


# test for u()
test_that('test that u computes calculations correctly', {
  
  tk <- c(0.000, 1.052)
  h_tk <- c(-3.751, -3.682)
  hprime_tk <- c(0.069, 0.062)
  x_1 <- 0.3
  
  # case when zk is an empty vector
  zk_1 <- c()
  
  # manual calculation
  u_man_1 <- -3.751 + (x_1 - 0.000) * 0.069
  
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


#Test for l() function
test_that('test the l function',{
  
  ## case when Tk only has one element.
  
  x_1 <- 0.3
  zk_1 <- c(0.562)
  tk_1 <- c(0.000)
  h_tk_1 <- c(-3.751,-3.682,-3.518)
  hprime_tk_1 <- c(0.069,0.062,0.053)
  
  l_man_1 <- -3.751 + (x_1-0.000)*0.069
  l_val_1 <- l(x_1,zk_1,tk_1,h_tk_1,hprime_tk_1)
  
  expect_equal(l_1,l_val_1)
  
  
  ## case when x < Tk[1] or x > Tk[k]
  
  
  x_2 <- -1
  zk_2 <- c(0.562)
  tk_2 <- c(0.000,1.052)
  h_tk_2 <- c(-3.751,-3.682,-3.518)
  hprime_tk_2 <- c(0.069,0.062,0.053)
  
  l_2 <- -Inf
  l_val_2 <- l(x_2,zk_2,tk_2,h_tk_2,hprime_tk_2)
  expect_equal(l_2,l_val_2)
  
  
  x_3 <- 1.6
  zk_3 <- c(0.562)
  tk_3 <- c(0.000,1.052)
  h_tk_3 <- c(-3.751,-3.682,-3.518)
  hprime_tk_3 <- c(0.069,0.062,0.053)
  
  l_3 <- -Inf
  l_val_3 <- l(x_3,zk_3,tk_3,h_tk_3,hprime_tk_3)
  expect_equal(l_3,l_val_3)
  
  ## case When x is on the left hand side of zk but bigger than tk[1]
  
  x_4 <- 0.4
  zk_4 <- c(0.562)
  tk_4 <- c(0.000,1.052)
  h_tk_4 <- c(-3.751,-3.682,-3.518)
  hprime_tk_4 <- c(0.069,0.062,0.053)
  
  l_4 <- ((1.052-x_4)*(-3.751) + (x_4-0.000)*(-3.682))/(1.052-0.000)
  l_val_4 <- l(x_4,zk_4,tk_4,h_tk_4,hprime_tk_4)
  expect_equal(l_4,l_val_4)
  
  ## case When x in [Tk_j,Tk_j+1]
  
  
  x_5 <- 1.82
  zk_5 <- c(0.562)
  tk_5 <- c(0.000,1.052,2.105)
  h_tk_5 <- c(-3.751,-3.682,-3.518)
  hprime_tk_5 <- c(0.069,0.062,0.053)
  
  l_5 <- ((2.105-x_5)*(-3.682) + (x_5-1.052)*(-3.518))/(2.105-1.052)
  l_val_5 <- l(x_5,zk,tk_5,h_tk_5,hprime_tk_5)
  expect_equal(l_5,l_val_5)
  
})


