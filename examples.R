# EXAMPLES - data simulation using NCC package

devtools::install_github("pavlakrotka/NCC", build_vignettes = TRUE)
library(NCC)

set.seed(123)

## Continuous endpoints - equal linear time trends; alternative hypothesis

cont_eq <- datasim_cont(n_total = 1000, num_arms = 3, d = 120, theta = rep(0.25, 3), lambda = rep(0.15, 4), sigma = 1, trend = "linear")

## Continuous endpoints - different linear time trends; alternative hypothesis

cont_diff <- datasim_cont(n_total = 1000, num_arms = 3, d = 120, theta = rep(0.25, 3), lambda = c(0.1, 0.15, 0.1, 0.1), sigma = 1, trend = "stepwise")

## Binary endpoints - equal linear time trends; alternative hypothesis

bin_eq <- datasim_bin(n_total = 1000, num_arms = 3, d = 120, p0 = 0.7, OR = rep(1.8, 3), lambda = rep(0.15, 4), trend = "linear")

## Binary endpoints - different linear time trends; alternative

bin_diff <- datasim_bin(n_total = 1000, num_arms = 3, d = 120, p0 = 0.7, OR = rep(1.8, 3), lambda = c(0.1, 0.15, 0.1, 0.1), trend = "linear")