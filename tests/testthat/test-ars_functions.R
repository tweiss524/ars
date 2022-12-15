
# test for log_concave()
test_that('log_concave errors for non-decreasing vector',{

  # create a non-log-concave vector, i.e., non-decreasing
  non_concave <- c(3, 1, 5, -1, 0)
  expect_error(check_log_concave(non_concave), "Function is not log-concave")

})


# test for calc_z()
test_that('test that calc_z computes calculation correctly',{

  # case when there is only 1 value in Tk
  tk_1 <- 1
  h_tk1 <- -3.751
  hprime_tk1 <- 0.069

  calc_z_with_1_value <- calc_z(tk_1, h_Tk1, hprime_Tk1)
  expect_equal(calc_z_with_1_value, c())

  # case where Tk has more than 1 value
  tk_2 <- c(0.000, 1.052)
  h_tk2 <- c(3.751, -3.682)
  hprime_tk2 <- c(0.069, 0.062)

  calc_z_normal_case <- calc_z(tk_2, h_tk2, hprime_tk2)

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
  u_1 <- -3.751 + (x_1 - 0.000) * 0.069
  u_val_1 <- u(x_1, zk_1, tk, h_tk, hprime_tk)

  expect_equal(u_1, u_val_1)


  # case when zk has more than 1 value
  zk_2 <- c(0.526, 0.777)

  # case when x is on the left hand side of zk
  u_2 <- -3.751 + (x_1-0.000)*0.069
  u_val_2 <- u(x_1, zk_2, tk, h_tk, hprime_tk)
  expect_equal(u_2, u_val_2)

  # case when x is on the right hand side of zk
  x_2 <- 0.6
  u_3 <- -3.682 + (x_2 - 1.052) * 0.062
  u_val_3 <- u(x_2, zk_2, tk, h_tk, hprime_tk)
  expect_equal(u_3, u_val_3)

})

# #Test for l() function
# test_that('test the l function',{
#
#
#   h_tk <- c(-3.751,-3.682,-3.518)
#   hprime_tk <- c(0.069,0.062,0.053)
#   x_1 <- 0.3
#
#   ## case when Tk only has one element.
#   tk_1 <- c(0.000)
#   zk <- c(0.562)
#   l_1 <- -3.751 + (x_1-0.000)*0.069
#   l_val_1 <- l(x_1,zk,tk_1,h_tk,hprime_tk)
#
#   expect_equal(l_1,l_val_1)
#
#
#   ## case when x < Tk[1] or x > Tk[k]
#   tk_2 <- c(0.000,1.052)
#
#   x_2 <- -1
#   l_2 <- -Inf
#   l_val_2 <- l(x_2,zk,tk_2,h_tk,hprime_tk)
#   expect_equal(l_2,l_val_2)
#
#
#   x_3 <- 1.6
#   l_3 <- -Inf
#   l_val_3 <- l(x_3,zk,tk_2,h_tk,hprime_tk)
#   expect_equal(l_3,l_val_3)
#
#   ## case When x is on the left hand side of zk but bigger than tk[1]
#   x_4 <- 0.4
#   l_4 <- ((1.052-x_4)*(-3.751) + (x_4-0.000)*(-3.682))/(1.052-0.000)
#   l_val_4 <- l(x_4,zk,tk_2,h_tk,hprime_tk)
#   expect_equal(l_4,l_val_4)
#
#   ## case When x is on the right hand side of zk but less than tk[3]
#
#   tk_3 <- c(0.000,1.052,2.105)
#   x_5 <- 0.6
#   l_5 <- ((1.679-x_5)*(-3.682) + (x_5-1.052)*(-3.518))/(1.679-1.052)
#   l_val_5 <- l(x_5,zk,tk_3,h_tk,hprime_tk)
#   expect_equal(l_5,l_val_5)
#
#
#
# })
