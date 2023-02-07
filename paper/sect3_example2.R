##################################################################################
# Code for the paper
# ``NCC: An R-package for analysis and simulation of platform trials with non-concurrent controls'' by Pavla Krotka et al. (2023)
# Example 2 - Sect 3
##################################################################################

library(NCC)
library(ggplot2)

lambda_values <- rep(seq(-0.15, 0.15, length.out = 9), 2)
sim_scenarios <- data.frame(num_arms = 4, 
                            n_arm = 250, 
                            d1 = 250*0,
                            d2 = 250*1,
                            d3 = 250*2,
                            d4 = 250*3,
                            period_blocks = 2, 
                            mu0 = 0,
                            sigma = 1,
                            theta1 = 0,
                            theta2 = 0,
                            theta3 = 0,
                            theta4 = 0,
                            lambda0 = rep(seq(-0.15, 0.15, length.out = 9), 2),
                            lambda1 = rep(seq(-0.15, 0.15, length.out = 9), 2),
                            lambda2 = rep(seq(-0.15, 0.15, length.out = 9), 2),
                            lambda3 = rep(seq(-0.15, 0.15, length.out = 9), 2),
                            lambda4 = rep(seq(-0.15, 0.15, length.out = 9), 2),
                            trend = c(rep("linear", 9), rep("stepwise_2", 9)),
                            alpha = 0.025,
                            ncc = TRUE)

head(sim_scenarios)

set.seed(1234)
sim_results <- sim_study_par(nsim = 1000, scenarios = sim_scenarios, arms = 4, 
                             models = c("fixmodel", "sepmodel", "poolmodel"), endpoint = "cont")

head(sim_results)

ggplot(sim_results, aes(x=lambda0, y=reject_h0, color=model)) +
  geom_point() +
  geom_line() +
  facet_grid(~ trend) +
  geom_hline(aes(yintercept = 0.025), linetype = "dotted") +
  labs(x="Strength of time trend", y="Type I error", color="Analysis approach") +
  theme_bw()

ggplot(sim_results, aes(x=lambda0, y=bias, color=model)) +
  geom_point() +
  geom_line() +
  facet_grid(~ trend) +
  geom_hline(aes(yintercept = 0), linetype = "dotted") +
  labs(x="Strength of time trend", y="Bias", color="Analysis approach") +
  theme_bw()

ggplot(sim_results, aes(x=lambda0, y=MSE, color=model)) +
  geom_point() +
  geom_line() +
  facet_grid(~ trend) +
  labs(x="Strength of time trend", y="MSE", color="Analysis approach") +
  theme_bw()