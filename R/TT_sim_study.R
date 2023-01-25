#' Wrapper function performing simulation studies (treatment-treatment comparisons) for a given set of scenarios
#'
#' @description This function performs a simulation study for a given set of scenarios, analyzing simulated data using the fixed effect model and indirect comparison. It compares the indicated treatment arms.
#'
#' @param nsim Number of replications
#' @param scenarios Data frame containing all parameters for scenarios that should be simulated
#' @param arms Indicator of the treatment arms to be compared (vector of length 2)
#' @param models Vector with models that should be used for the analysis. Default=c("fixmodel", "indirect")
#' @param endpoint Endpoint indicator. "cont" for continuous endpoints, "bin" for binary endpoints
#'
#' @keywords internal
#'
#' @export
#'
#'
#' @return Data frame with all considered scenarios and corresponding results
#' @author Pavla Krotka


TT_sim_study <- function(nsim, scenarios, arms, models = c("fixmodel", "indirect"), endpoint){

  models <- sort(models)

  if(endpoint=="cont"){

    result <- NULL

    num_models <- length(models)

    for(i in 1:dim(scenarios)[1]){

      d_i <- as.numeric(scenarios[i, grepl("^d\\d", names(scenarios))])
      theta_i <- as.numeric(scenarios[i, grepl("^theta\\d", names(scenarios))])
      lambda_i <- as.numeric(scenarios[i, grepl("^lambda\\d", names(scenarios))])

      db <- replicate(nsim,
                      TT_all_models(data = datasim_cont(n_arm = scenarios$n_arm[i],
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
                                    arms = arms,
                                    models = models,
                                    endpoint = endpoint,
                                    alpha = scenarios$alpha[i]))

      result_i <- cbind(scenarios[i,],
                        compared_arms = paste0(arms[1], " vs. ", arms[2]),
                        model = models,
                        reject_h0 = rowMeans(matrix(as.logical(unlist(unname(db))), ncol = nsim), na.rm = TRUE)) # get power/T1E

      result <- rbind(result, result_i)
    }
  }


  if(endpoint=="bin"){

    result <- NULL

    num_models <- length(models)

    for(i in 1:dim(scenarios)[1]){

      d_i <- as.numeric(scenarios[i, grepl("^d\\d", names(scenarios))])
      OR_i <- as.numeric(scenarios[i, grepl("^OR\\d", names(scenarios))])
      lambda_i <- as.numeric(scenarios[i, grepl("^lambda\\d", names(scenarios))])

      db <- replicate(nsim,
                      TT_all_models(data = datasim_bin(n_arm = scenarios$n_arm[i],
                                                       num_arms = scenarios$num_arms[i],
                                                       d = d_i,
                                                       period_blocks = scenarios$period_blocks[i],
                                                       p0 = scenarios$p0[i],
                                                       OR = OR_i,
                                                       lambda = lambda_i,
                                                       trend = scenarios$trend[i],
                                                       N_peak = scenarios$N_peak[i],
                                                       full = FALSE),
                                    arms = arms,
                                    models = models,
                                    endpoint = endpoint,
                                    alpha = scenarios$alpha[i]))

      result_i <- cbind(scenarios[i,],
                        compared_arms = paste0(arms[1], " vs. ", arms[2]),
                        model = models,
                        reject_h0 = rowMeans(matrix(as.logical(unlist(unname(db))), ncol = nsim), na.rm = TRUE)) # get power/T1E

      result <- rbind(result, result_i)

    }
  }
  return(result)
}



