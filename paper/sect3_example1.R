##################################################################################
# Code for the paper
# ``NCC: An R-package for analysis and simulation of platform trials with non-concurrent controls'' by Pavla Krotka et al. (2023)
# Example - Sect 3
##################################################################################

library(NCC)

set.seed(2346)
data <- datasim_bin(num_arms = 3, n_arm = 100, d = c(0, 100, 250),
                    p0 = 0.7, OR = rep(1.8, 3), lambda = rep(0.15, 4), trend="stepwise")

head(data)
plot_trial(data$treatment)

fmodel <- fixmodel_bin(data=data,arm=3)
summary(fmodel$model)

MAPprior_bin(data=data,arm=3)

timemachine_bin(data=data,arm=3)
