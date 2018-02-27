context("test-pm_plot.R")




test_that("pm_plot works with variables in data frame and parent env", {
  a <- rep(c(1,1,1,0),100)
  b <- rep(c(2,0,3,0),100)

  df <- data.frame(a,a,a,b,a,b,a,b)


  expect_is(pm_plot(df), "ggplot")
  expect_is(pm_plot(df,show_cluster = 3), "ggplot")
  expect_is(pm_plot(df,plot_type = "bar"), "ggplot")
  expect_is(pm_plot(df,show_cluster = 3,plot_type = "bar"), "ggplot")
  expect_is(pm_plot(df,use_log = FALSE), "ggplot")
})
