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
#' @importFrom parallel makeCluster
#' @importFrom parallel stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach
#' @importFrom foreach %dopar%
#' @importFrom iterators icount
#'
#' @keywords internal
#'
#' @export
#'
#'
#' @return Data frame with all considered scenarios and corresponding results
#' @author Pavla Krotka


TT_sim_study_par <- function(nsim, scenarios, arms, models = c("fixmodel", "indirect"), endpoint, perc_cores=0.9){

  cores <- detectCores()
  cl <- makeCluster(floor(cores[1]*perc_cores)) # not to overload your computer
  registerDoParallel(cl)

  models <- sort(models)

  if(endpoint=="cont"){

    result <- NULL

    num_models <- length(models)

    for(i in 1:dim(scenarios)[1]){

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
                                 alpha = scenarios$alpha[i])

                    }



      result_i <- cbind(scenarios[i,],
                        compared_arms = paste0(arms[1], " vs. ", arms[2]),
                        model = models,
                        reject_h0 = rowMeans(matrix(as.logical(unlist(unname(db))), ncol = nsim), na.rm = TRUE)) # get power/T1E

      result <- rbind(result, result_i)
    }

    stopCluster(cl)
    gc()

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
                                 alpha = scenarios$alpha[i])
                    }



      result_i <- cbind(scenarios[i,],
                        compared_arms = paste0(arms[1], " vs. ", arms[2]),
                        model = models,
                        reject_h0 = rowMeans(matrix(as.logical(unlist(unname(db))), ncol = nsim), na.rm = TRUE)) # get power/T1E

      result <- rbind(result, result_i)

    }

    stopCluster(cl)
    gc()

  }
  return(result)
}



