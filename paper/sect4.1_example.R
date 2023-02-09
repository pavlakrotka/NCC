##################################################################################
# Code for the paper
# "NCC: An R-package for analysis and simulation of platform trials with non-concurrent controls" by Pavla Krotka et al. (2023)
# Example - Section 4
##################################################################################

library(NCC)

set.seed(123)
trial_data <- datasim_bin(num_arms = 3, n_arm = 100, d = c(0, 100, 250),
                    p0 = 0.7, OR = rep(1.8, 3), lambda = rep(0.15, 4), trend="stepwise")

head(trial_data)

plot_trial(trial_data$treatment)
#ggplot2::ggsave("plot_trial.png", width = 7, height = 4)

fixmodel_bin(data=trial_data, arm=3, alpha=0.025)
summary(fixmodel_bin(data=trial_data, arm=3)$model)

MAPprior_bin(data=trial_data, arm=3, alpha=0.025)

timemachine_bin(data=trial_data, arm=3, alpha=0.025)
