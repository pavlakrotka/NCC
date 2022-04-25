# EXAMPLES - data simulation using NCC package

devtools::install_github("pavlakrotka/NCC", build_vignettes = TRUE)
library(NCC)

set.seed(123)

## Continuous endpoints - equal linear time trends; alternative hypothesis

cont_eq <- datasim_cont(num_arms = 3, n_arm = 100, d = c(0, 100, 250), theta = rep(0.25, 3), lambda = rep(0.15, 4), sigma = 1, trend = "linear")

## Continuous endpoints - different linear time trends; alternative hypothesis

cont_diff <- datasim_cont(num_arms = 3, n_arm = 100, d = c(0, 100, 250), theta = rep(0.25, 3), lambda = c(0.1, 0.15, 0.1, 0.1), sigma = 1, trend = "stepwise")

## Binary endpoints - equal linear time trends; alternative hypothesis

bin_eq <- datasim_bin(num_arms = 3, n_arm = 100, d = c(0, 100, 250), p0 = 0.7, OR = rep(1.8, 3), lambda = rep(0.15, 4), trend = "linear")

## Binary endpoints - different linear time trends; alternative hypothesis

bin_diff <- datasim_bin(num_arms = 3, n_arm = 100, d = c(0, 100, 250), p0 = 0.7, OR = rep(1.8, 3), lambda = c(0.1, 0.15, 0.1, 0.1), trend = "linear")
