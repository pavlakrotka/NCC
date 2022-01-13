
# # Copy in R folder the functions of the r package 
# setwd("C:/Users/mbofi/Nextcloud/GitKraken/NCC")
# devtools::document()
# devtools::load_all()
# 
# # Build & check the package
# devtools::build(pkg = "C:/Users/mbofi/Nextcloud/GitKraken/NCC", path = NULL, binary = FALSE, manual = TRUE)
# devtools::check_built(path = "C:/Users/mbofi/Nextcloud/GitKraken/NCC", cran=TRUE, manual = TRUE)
# devtools::build_manual(pkg = "C:/Users/mbofi/Nextcloud/GitKraken/NCC", path = "C:/Users/mbofi/Nextcloud/GitKraken/NCC")


devtools::install_github("pavlakrotka/NCC",auth_token = "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", force=TRUE)

library(NCC)

?data_sim

# Example

A0 <- c(10, 10, 10)
A1 <- c(20, 10, NA)
A2 <- c(NA, 10, 20)

SS_matrix <- matrix(c(A0, A1, A2), nrow = 3, byrow = T)
SS_matrix

alloc_ratios <- matrix(c(1,1,1,
                         2,1,0,
                         0,1,2), ncol = 3, byrow = T)
alloc_ratios

test <- data_sim(SS_matrix = SS_matrix, block_sizes = c(6,6,6), alloc_ratios = alloc_ratios,
                 mu0 = 0, delta = c(0.25,0.25), p0, OR, lambda = c(0,0,0), sigma = 1, N_peak, trend = "linear", trend_param, endpoint = "continuous")

head(test)
