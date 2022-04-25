#' Wrapper function performing simulation studies for a given set of scenarios
#'
#' @description Performs a simulation study for a given set of scenarios, analyzing simulated data using the fixed effect model, pooled and separate analyses, and timemachine and MAP prior approach. Performs inference for all treatment arms in the trial except for the first one
#'
#' @param nsim Number of replications
#' @param scenarios Data frame containing all parameters for scenarios that should be simulated
#' @param endpoint Endpoint indicator. "cont" for continuous endpoints, "bin" for binary endpoints
#'
#' @keywords internal
#'
#' @export
#'
#'
#' @return List containing an indicator whether the null hypothesis was rejected or not for the investigated treatment for all models
#' @author Pavla Krotka


sim_study <- function(nsim, scenarios, endpoint){

  if(endpoint=="cont"){

    result <- NULL

    num_models <- 4

    for(i in 1:dim(scenarios)[1]){

      d_i <- as.numeric(scenarios[i, grepl("^d\\d", names(scenarios))])
      theta_i <- as.numeric(scenarios[i, grepl("^theta\\d", names(scenarios))])
      lambda_i <- as.numeric(scenarios[i, grepl("^lambda\\d", names(scenarios))])

      db <- replicate(nsim,
                      all_models(data = datasim_cont(n_arm = scenarios$n_arm[i],
                                                     num_arms = scenarios$num_arms[i],
                                                     d = d_i,
                                                     period_blocks = scenarios$period_blocks[i],
                                                     mu0 = scenarios$mu0[i],
                                                     theta = theta_i,
                                                     lambda = lambda_i,
                                                     sigma = scenarios$sigma[i],
                                                     trend = scenarios$trend[i],
                                                     N_peak = scenarios$N_peak[i],
                                                     full = FALSE),
                                 alpha = scenarios$alpha[i]))

      result_i <- cbind(scenarios[i,],
                        study_arm = rep(c(2:((nrow(db)/num_models)+1)), each = num_models), # inference for all arms starting with 2
                        model = c("fix", "sep", "pool", "timemachine"),
                        reject_h0 = rowMeans(matrix(unlist(unname(db)), ncol = nsim))) # get power/T1E

      result <- rbind(result, result_i)
    }
  }


  if(endpoint=="bin"){

    result <- NULL

    num_models <- 5

    for(i in 1:dim(scenarios)[1]){

      d_i <- as.numeric(scenarios[i, grepl("^d\\d", names(scenarios))])
      OR_i <- as.numeric(scenarios[i, grepl("^OR\\d", names(scenarios))])
      lambda_i <- as.numeric(scenarios[i, grepl("^lambda\\d", names(scenarios))])

      db <- replicate(nsim,
                      all_models(data = datasim_bin(n_arm = scenarios$n_arm[i],
                                                    num_arms = scenarios$num_arms[i],
                                                    d = d_i,
                                                    period_blocks = scenarios$period_blocks[i],
                                                    p0 = scenarios$p0[i],
                                                    OR = OR_i,
                                                    lambda = lambda_i,
                                                    trend = scenarios$trend[i],
                                                    N_peak = scenarios$N_peak[i],
                                                    full = FALSE),
                                 alpha = scenarios$alpha[i]))

      result_i <- cbind(scenarios[i,],
                        study_arm = rep(c(2:((nrow(db)/num_models)+1)), each = num_models), # inference for all arms starting with 2
                        model = c("fix", "sep", "pool", "timemachine", "MAPprior"),
                        reject_h0 = rowMeans(matrix(unlist(unname(db)), ncol = nsim))) # get power/T1E

      result <- rbind(result, result_i)

    }
  }
  return(result)
}



