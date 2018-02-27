library(testthat)
library(pasaR)

test_check("pasaR")
use_coverage(pkg = ".", type = c("codecov"))
