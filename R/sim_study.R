#' Wrapper function performing simulation studies for a given set of scenarios
#'
#' @description Performs a simulation study for a given set of scenarios, analyzing simulated data using the fixed effect model, pooled and separate analyses, and timemachine and MAP prior approach. Performs inference for all treatment arms in the trial except for the first one
#'
#' @param nsim Number of replications
#' @param scenarios Data frame containing all parameters for scenarios that should be simulated
#' @param arms Vector with treatment arms to perform inference on. Default - all arms except the first one
#' @param models Vector with models that should be used for the analysis. Default=c("fixmodel", "sepmodel", "poolmodel", "timemachine", "MAPprior")
#' @param endpoint Endpoint indicator. "cont" for continuous endpoints, "bin" for binary endpoints
#'
#' @keywords internal
#'
#' @export
#'
#'
#' @return Data frame with all considered scenarios and corresponding results
#' @author Pavla Krotka


sim_study <- function(nsim, scenarios, arms, models = c("fixmodel", "sepmodel", "poolmodel", "timemachine", "mixmodel", "MAP_rjags"), endpoint){

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
                                 arms = arms,
                                 models = models,
                                 endpoint = endpoint,
                                 alpha = scenarios$alpha[i],
                                 unit_size = scenarios$unit_size[i],
                                 ncc = scenarios$ncc[i],
                                 opt = scenarios$opt[i],
                                 prior_prec_tau = scenarios$prior_prec_tau[i],
                                 n.samples = scenarios$n.samples[i],
                                 n.chains = scenarios$n.chains[i],
                                 n.iter = scenarios$n.iter[i],
                                 n.adapt = scenarios$n.adapt[i],
                                 robustify = scenarios$robustify[i],
                                 weight = scenarios$weight[i],
                                 ci = scenarios$ci[i],
                                 prec_delta = scenarios$prec_delta[i],
                                 prec_gamma = scenarios$prec_gamma[i],
                                 tau_a = scenarios$tau_a[i],
                                 tau_b = scenarios$tau_b[i],
                                 prec_a = scenarios$prec_a[i],
                                 prec_b = scenarios$prec_b[i],
                                 bucket_size = scenarios$bucket_size[i],
                                 smoothing_basis = scenarios$smoothing_basis[i],
                                 basis_dim = scenarios$basis_dim[i],
                                 gam_method = scenarios$gam_method[i],
                                 bs_degree = scenarios$bs_degree[i]))

      result_i <- cbind(scenarios[i,],
                        study_arm = rep(arms, each = num_models),
                        model = models,
                        reject_h0 = rowMeans(matrix(as.logical(unlist(unname(db))), ncol = nsim), na.rm = TRUE), # get power/T1E
                        failed = rowSums(is.na(matrix(as.logical(unlist(unname(db))), ncol = nsim))),
                        nsim = nsim)

      result <- rbind(result, result_i)

      print(paste0("Scenario ", i, "/", dim(scenarios)[1], " done"))
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
                                 arms = arms,
                                 models = models,
                                 endpoint = endpoint,
                                 alpha = scenarios$alpha[i],
                                 unit_size = scenarios$unit_size[i],
                                 ncc = scenarios$ncc[i],
                                 opt = scenarios$opt[i],
                                 prior_prec_tau = scenarios$prior_prec_tau[i],
                                 n.samples = scenarios$n.samples[i],
                                 n.chains = scenarios$n.chains[i],
                                 n.iter = scenarios$n.iter[i],
                                 n.adapt = scenarios$n.adapt[i],
                                 robustify = scenarios$robustify[i],
                                 weight = scenarios$weight[i],
                                 ci = scenarios$ci[i],
                                 prec_delta = scenarios$prec_delta[i],
                                 prec_gamma = scenarios$prec_gamma[i],
                                 tau_a = scenarios$tau_a[i],
                                 tau_b = scenarios$tau_b[i],
                                 prec_a = scenarios$prec_a[i],
                                 prec_b = scenarios$prec_b[i],
                                 bucket_size = scenarios$bucket_size[i],
                                 smoothing_basis = scenarios$smoothing_basis[i],
                                 basis_dim = scenarios$basis_dim[i],
                                 gam_method = scenarios$gam_method[i],
                                 bs_degree = scenarios$bs_degree[i]))

      result_i <- cbind(scenarios[i,],
                        study_arm = rep(arms, each = num_models),
                        model = models,
                        reject_h0 = rowMeans(matrix(as.logical(unlist(unname(db))), ncol = nsim), na.rm = TRUE), # get power/T1E
                        failed = rowSums(is.na(matrix(as.logical(unlist(unname(db))), ncol = nsim))),
                        nsim = nsim)

      result <- rbind(result, result_i)

      print(paste0("Scenario ", i, "/", dim(scenarios)[1], " done"))

    }
  }
  return(result)
}



