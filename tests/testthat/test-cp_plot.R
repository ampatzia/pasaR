context("test-cp_plot.R")

test_that("cp_plot works with variables in data frame and parent env", {
  a <- rep(c(1,1,1,0),100)
  b <- rep(c(2,0,3,0),100)

  df <- data.frame(a,a,a,b,a,b,a,b)


  expect_is(cp_plot(df), "ggplot")
  expect_is(cp_plot(df,show_cluster = 3), "ggplot")
  expect_is(cp_plot(df,plot_type = "bar"), "ggplot")
  expect_is(cp_plot(df,show_cluster = 3,plot_type = "bar"), "ggplot")
  expect_is(cp_plot(df,use_log = FALSE), "ggplot")
})
