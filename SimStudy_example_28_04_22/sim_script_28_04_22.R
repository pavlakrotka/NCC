#devtools::install_github("pavlakrotka/NCC", build = TRUE, force=T)
library(NCC)
library(tidyverse)

n_sim <- 1000
set.seed(123)

# Scenario II

#get_ss_matrix(num_arms = 3, n_arm = 100, d = c(0,0,0))
#get_ss_matrix(num_arms = 3, n_arm = 100, d = c(0,100,200))
#get_ss_matrix(num_arms = 3, n_arm = 100, d = c(0,200,400))


scenario_ii_pow <- data.frame(num_arms = 3, 
                              n_arm = 100, 
                              d1 = 0,
                              d2 = c(0, 100, 200),
                              d3 = c(0, 200, 400),
                              period_blocks = 2, 
                              p0 = 0.3, 
                              OR1 = 1.8,
                              OR2 = 1.8,
                              OR3 = 1.8,
                              lambda0 = 0.3, 
                              lambda1 = c(rep(0.3, 3), rep(0.5, 3)),
                              lambda2 = 0.3,
                              lambda3 = 0.3,
                              trend = "stepwise",
                              alpha = 0.1) %>%
  mutate(timetrend = ifelse(lambda0==lambda1 & lambda1==lambda2 & lambda2==lambda3, "EQ", "DIFF"),
         hypothesis = ifelse(OR2==1 & OR3==1, "H0", "H1"))

results_ii_pow <- sim_study(nsim = n_sim, scenarios = scenario_ii_pow, endpoint = "bin")
write_csv(results_ii_pow, "results/results_ii_pow.csv")




scenario_ii_alpha <- data.frame(num_arms = 3, 
                                n_arm = 100, 
                                d1 = 0,
                                d2 = c(0, 100, 200),
                                d3 = c(0, 200, 400),
                                period_blocks = 2, 
                                p0 = 0.3, 
                                OR1 = 1.8,
                                OR2 = 1,
                                OR3 = 1,
                                lambda0 = 0.3, 
                                lambda1 = c(rep(0.3, 3), rep(0.5, 3)),
                                lambda2 = 0.3,
                                lambda3 = 0.3,
                                trend = "stepwise",
                                alpha = 0.1) %>%
  mutate(timetrend = ifelse(lambda0==lambda1 & lambda1==lambda2 & lambda2==lambda3, "EQ", "DIFF"),
         hypothesis = ifelse(OR2==1 & OR3==1, "H0", "H1"))

results_ii_alpha <- sim_study(nsim = n_sim, scenarios = scenario_ii_alpha, endpoint = "bin")
write_csv(results_ii_alpha, "results/results_ii_alpha.csv")

# Scenario III

## get_ss_matrix(num_arms = 4, n_arm = 100, d = c(0,100,100,200))

#get_ss_matrix(num_arms = 4, n_arm = 100, d = c(0,100,150,200))

#get_ss_matrix(num_arms = 4, n_arm = 100, d = c(0,100,200,200))



scenario_iii_pow <- data.frame(num_arms = 4, 
                               n_arm = 100, 
                               d1 = 0,
                               d2 = 100,
                               d3 = c(150, 200),
                               d4 = 200,
                               period_blocks = 2, 
                               p0 = 0.3, 
                               OR1 = 1.8,
                               OR2 = 1.8,
                               OR3 = 1.8,
                               OR4 = 1.8,
                               lambda0 = 0.3, 
                               lambda1 = c(rep(0.3, 2), rep(0.5, 2)),
                               lambda2 = 0.3,
                               lambda3 = 0.3,
                               lambda4 = 0.3,
                               trend = "stepwise",
                               alpha = 0.1) %>%
  mutate(timetrend = ifelse(lambda0==lambda1 & lambda1==lambda2 & lambda2==lambda3 & lambda3==lambda4, "EQ", "DIFF"),
         hypothesis = ifelse(OR2==1 & OR3==1 & OR4==1, "H0", "H1"))

results_iii_pow <- sim_study(nsim = n_sim, scenarios = scenario_iii_pow, endpoint = "bin")
write_csv(results_iii_pow, "results/results_iii_pow.csv")



scenario_iii_alpha <- data.frame(num_arms = 4, 
                                 n_arm = 100, 
                                 d1 = 0,
                                 d2 = 100,
                                 d3 = c(150, 200),
                                 d4 = 200,
                                 period_blocks = 2, 
                                 p0 = 0.3, 
                                 OR1 = 1.8,
                                 OR2 = 1,
                                 OR3 = 1,
                                 OR4 = 1,
                                 lambda0 = 0.3, 
                                 lambda1 = c(rep(0.3, 2), rep(0.5, 2)),
                                 lambda2 = 0.3,
                                 lambda3 = 0.3,
                                 lambda4 = 0.3,
                                 trend = "stepwise",
                                 alpha = 0.1) %>%
  mutate(timetrend = ifelse(lambda0==lambda1 & lambda1==lambda2 & lambda2==lambda3 & lambda3==lambda4, "EQ", "DIFF"),
         hypothesis = ifelse(OR2==1 & OR3==1 & OR4==1, "H0", "H1"))

results_iii_alpha <- sim_study(nsim = n_sim, scenarios = scenario_iii_alpha, endpoint = "bin")
write_csv(results_iii_alpha, "results/results_iii_alpha.csv")

