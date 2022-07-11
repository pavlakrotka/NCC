#' Wrapper function performing simulation studies for a given set of scenarios (parallelized)
#'
#' @description Performs a simulation study for a given set of scenarios, analyzing simulated data using the fixed effect model, pooled and separate analyses, and timemachine and MAP prior approach. Performs inference for all treatment arms in the trial except for the first one
#'
#' @param nsim Number of replications
#' @param scenarios Data frame containing all parameters for scenarios that should be simulated
#' @param arms Vector with treatment arms to perform inference on. Default - all arms except the first one
#' @param models Vector with models that should be used for the analysis. Default=c("fixmodel", "sepmodel", "poolmodel", "timemachine", "MAPprior")
#' @param endpoint Endpoint indicator. "cont" for continuous endpoints, "bin" for binary endpoints
#'
#' @importFrom parallel detectCores
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach
#' @importFrom foreach %dopar%
#'
#' @keywords internal
#'
#' @export
#'
#'
#' @return Data frame with all considered scenarios and corresponding results
#' @author Pavla Krotka


sim_study_par <- function(nsim, scenarios, arms, models = c("fixmodel", "sepmodel", "poolmodel", "timemachine", "mixmodel", "MAPprior"), endpoint){

  cores <- detectCores()-1 # not to overload your computer
  registerDoParallel(cores)

  if (endpoint=="cont") {
    models <- models[models!="MAPprior"] # not implemented yet
  }

  if (endpoint=="bin") {
    models <- models[models!="mixmodel"] # not implemented yet
  }

  models <- sort(models)

  if(endpoint=="cont"){

    result <- NULL

    num_models <- length(models)

    for(i in 1:dim(scenarios)[1]){

      if(missing(arms)){
        arms <- c(2:scenarios[i,]$num_arms)
      }

      arms <- sort(arms)

      d_i <- as.numeric(scenarios[i, grepl("^d\\d", names(scenarios))])
      theta_i <- as.numeric(scenarios[i, grepl("^theta\\d", names(scenarios))])
      lambda_i <- as.numeric(scenarios[i, grepl("^lambda\\d", names(scenarios))])

      # db <- replicate(nsim,
      #                 all_models(data = datasim_cont(n_arm = scenarios$n_arm[i],
      #                                                num_arms = scenarios$num_arms[i],
      #                                                d = d_i,
      #                                                period_blocks = scenarios$period_blocks[i],
      #                                                mu0 = scenarios$mu0[i],
      #                                                theta = theta_i,
      #                                                lambda = lambda_i,
      #                                                sigma = scenarios$sigma[i],
      #                                                trend = scenarios$trend[i],
      #                                                N_peak = scenarios$N_peak[i],
      #                                                full = FALSE),
      #                            arms = arms,
      #                            models = models,
      #                            endpoint = endpoint,
      #                            alpha = scenarios$alpha[i]))

      db <- foreach(icount(nsim), .combine = cbind,
                    .packages = c("NCC")) %dopar% {
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
                                 arms = arms,
                                 models = models,
                                 endpoint = endpoint,
                                 alpha = scenarios$alpha[i])

                    }

      result_i <- cbind(scenarios[i,],
                        study_arm = rep(arms, each = num_models),
                        model = models,
                        reject_h0 = rowMeans(matrix(unlist(unname(db)), ncol = nsim))) # get power/T1E

      result <- rbind(result, result_i)
    }
  }


  if(endpoint=="bin"){

    result <- NULL

    num_models <- length(models)

    for(i in 1:dim(scenarios)[1]){

      if(missing(arms)){
        arms <- c(2:scenarios[i,]$num_arms)
      }

      arms <- sort(arms)

      d_i <- as.numeric(scenarios[i, grepl("^d\\d", names(scenarios))])
      OR_i <- as.numeric(scenarios[i, grepl("^OR\\d", names(scenarios))])
      lambda_i <- as.numeric(scenarios[i, grepl("^lambda\\d", names(scenarios))])

      # db <- replicate(nsim,
      #                 all_models(data = datasim_bin(n_arm = scenarios$n_arm[i],
      #                                               num_arms = scenarios$num_arms[i],
      #                                               d = d_i,
      #                                               period_blocks = scenarios$period_blocks[i],
      #                                               p0 = scenarios$p0[i],
      #                                               OR = OR_i,
      #                                               lambda = lambda_i,
      #                                               trend = scenarios$trend[i],
      #                                               N_peak = scenarios$N_peak[i],
      #                                               full = FALSE),
      #                            arms = arms,
      #                            models = models,
      #                            endpoint = endpoint,
      #                            alpha = scenarios$alpha[i]))


      db <- foreach(icount(nsim), .combine = cbind,
                    .packages = c("NCC")) %dopar% {
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
                                 arms = arms,
                                 models = models,
                                 endpoint = endpoint,
                                 alpha = scenarios$alpha[i])
                    }

      result_i <- cbind(scenarios[i,],
                        study_arm = rep(arms, each = num_models),
                        model = models,
                        reject_h0 = rowMeans(matrix(unlist(unname(db)), ncol = nsim))) # get power/T1E

      result <- rbind(result, result_i)

    }
  }
  return(result)
}



