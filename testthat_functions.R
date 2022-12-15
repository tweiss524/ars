
# Test for function log_concave
test_that('test wheather the log_concave function works',{
  
  # Creat a non-log-concave vector.
  non_concave <- c(3,1,5,-1,0)
  expect_error(check_log_concave(non_concave), "Function is not log-concave")
  
})


# Test for cal_z function.
test_that('test the cal_z function',{
  
  tk_1 <- c(1)
  tk_2 <- c(0.000,1.052)
  h_tk <- c(-3.751,-3.682)
  hprime_tk <- c(0.069,0.062)
  
  # case when there is only 1 value in Tk
  calc_z_with_1_value <- calc_z(tk_1,h_Tk,hprime_Tk)
  expect_equal(calc_z_with_1_value,c())
  
  # case when Tk that have more than 1 value.
  calc_z_normal_case <- calc_z(tk_2,h_tk,hprime_tk)
  z_val <- (-3.682-(-3.751)-1.052*0.062+0.000*0.069)/(0.069-0.062)
  expect_equal(calc_z_normal_case,z_val)
  
})


#Test for u function
test_that('test the u function',{
  
  tk <- c(0.000,1.052)
  h_tk <- c(-3.751,-3.682)
  hprime_tk <- c(0.069,0.062)
  x_1 <- 0.3
  
  ## case when zk is an empty vector.
  zk_1 <- c()
  u_1 <- -3.751 + (x_1-0.000)*0.069
  u_val_1 <- u(x_1,zk_1,tk,h_tk,hprime_tk)
  
  expect_equal(u_1,u_val_1)
  
  
  # case when Tk that have more than 1 value.
  zk_2 <- c(0.526)
  
  ## case When x is on the left hand side of zk
  u_2 <- -3.751 + (x_1-0.000)*0.069
  u_val_2 <- u(x_1,zk_2,tk,h_tk,hprime_tk)
  expect_equal(u_2,u_val_2)
  
  ## case When x is on the left hand side of zk
  x_2 <- 0.6
  u_3 <- -3.682 + (x_2-1.052)*0.062
  u_val_3 <- u(x_2,zk_2,tk,h_tk,hprime_tk)
  expect_equal(u_3,u_val_3)
  
})