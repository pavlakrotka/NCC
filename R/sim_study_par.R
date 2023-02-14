#' Wrapper function performing simulation studies for a given set of scenarios (parallelized on replication level)
#'
#' @description This function performs a simulation study for a given set of scenarios, analyzing simulated data using different models as indicated by the user. Performs inference for indicated experimental treatment arms. Simulates the probability to reject \eqn{H_0}, and the bias, as well as the mean squared error (MSE) of the treatment effect estimates based on a given number of replications.
#'
#' @param nsim Number of replications.
#' @param scenarios Data frame containing all parameters for scenarios that should be simulated.
#' @param arms Vector with treatment arms to perform inference on. These arms are compared to the control group. Default - all arms except the first one.
#' @param models Vector with models that should be used for the analysis. Default=c("fixmodel", "sepmodel", "poolmodel"). Available models for continuous endpoints are: 'fixmodel', 'fixmodel_cal', 'gam', 'MAPprior', 'mixmodel', 'mixmodel_cal', 'mixmodel_AR1', 'mixmodel_AR1_cal', 'piecewise', 'piecewise_cal', 'poolmodel', 'sepmodel', 'sepmodel_adj', 'splines', 'splines_cal', 'timemachine'. Available models for binary endpoints are: 'fixmodel', 'fixmodel_cal', 'MAPprior', 'poolmodel', 'sepmodel', 'sepmodel_adj', 'timemachine'.
#' @param endpoint Endpoint indicator. "cont" for continuous endpoints, "bin" for binary endpoints.
#' @param perc_cores What percentage of available cores should be used for the simulations. Default=0.9.
#'
#' @importFrom parallel detectCores
#' @importFrom parallelly availableCores
#' @importFrom parallel makeCluster
#' @importFrom parallel stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach
#' @importFrom foreach %dopar%
#' @importFrom iterators icount
#'
#' @export
#'
#' @examples
#'
#' # Create data frame with all parameters:
#' sim_scenarios <- data.frame(num_arms = 4,
#' n_arm = 250,
#' d1 = 250*0,
#' d2 = 250*1,
#' d3 = 250*2,
#' d4 = 250*3,
#' period_blocks = 2,
#' mu0 = 0,
#' sigma = 1,
#' theta1 = 0,
#' theta2 = 0,
#' theta3 = 0,
#' theta4 = 0,
#' lambda0 = rep(seq(-0.15, 0.15, length.out = 9), 2),
#' lambda1 = rep(seq(-0.15, 0.15, length.out = 9), 2),
#' lambda2 = rep(seq(-0.15, 0.15, length.out = 9), 2),
#' lambda3 = rep(seq(-0.15, 0.15, length.out = 9), 2),
#' lambda4 = rep(seq(-0.15, 0.15, length.out = 9), 2),
#' trend = c(rep("linear", 9), rep("stepwise_2", 9)),
#' alpha = 0.025,
#' ncc = TRUE)
#'
#' # Run simulation study:
#' # sim_results <- sim_study_par(nsim = 100, scenarios = sim_scenarios, arms = c(3, 4),
#' # models = c("fixmodel", "sepmodel", "poolmodel"), endpoint = "cont")
#'
#'
#' @return Data frame with all considered scenarios and corresponding results - the probability to reject \eqn{H_0}, and the bias, as well as the mean squared error (MSE) of the treatment effect estimates.
#' @author Pavla Krotka


sim_study_par <- function(nsim, scenarios, arms, models = c("fixmodel", "sepmodel", "poolmodel"), endpoint, perc_cores=0.9){

  if(!is.numeric(nsim) | length(nsim)!=1 | nsim<=1){
    stop("Number of replications (`nsim`) must be one number and must be larger than 1!")
  }

  if(!is.numeric(arms)){
    stop("Experimental treatment arms to be eveluated (`arms`) must be a numeric vector!")
  }

  if((endpoint %in% c("cont", "bin")==FALSE) | length(endpoint)!=1){
    stop("Endpoint indicator (`endpoint`) must be one of the following strings: 'cont', 'bin'!")
  }

  if(endpoint=="cont" & sum(models %in% c("fixmodel", "fixmodel_cal", "gam", "MAPprior",
                                          "mixmodel", "mixmodel_cal", "mixmodel_AR1", "mixmodel_AR1_cal",
                                          "piecewise", "piecewise_cal", "poolmodel", "sepmodel", "sepmodel_adj",
                                          "splines", "splines_cal", "timemachine")==FALSE)>0){
    stop("For continuous endpoints, only the following models are implemented: 'fixmodel', 'fixmodel_cal', 'gam', 'MAPprior',
                                          'mixmodel', 'mixmodel_cal', 'mixmodel_AR1', 'mixmodel_AR1_cal',
                                          'piecewise', 'piecewise_cal', 'poolmodel', 'sepmodel', 'sepmodel_adj',
                                          'splines', 'splines_cal', 'timemachine'.
         The argument `models` must contain only these strings!")
  }

  if(endpoint=="bin" & sum(models %in% c("fixmodel", "fixmodel_cal", "MAPprior", "poolmodel", "sepmodel", "sepmodel_adj", "timemachine")==FALSE)>0){
    stop("For binary endpoints, only the following models are implemented: 'fixmodel', 'fixmodel_cal', 'MAPprior', 'poolmodel', 'sepmodel', 'sepmodel_adj', 'timemachine'.
         The argument `models` must contain only these strings!")
  }

  if(!is.numeric(perc_cores) | length(perc_cores)!=1 | perc_cores>=1 | perc_cores<=0){
    stop("Percentage of cores to be used for simulations (`perc_cores`) must be one number between 0 and 1!")
  }






  print(paste0("Starting the simulations. ", dim(scenarios)[1], " scenarios will be simulated. Starting time: ", Sys.time()))

  cores <- availableCores()
  cl <- makeCluster(floor(unname(cores)*perc_cores)) # not to overload your computer
  registerDoParallel(cl)


  models <- sort(models)

  if(endpoint=="cont"){

    result <- NULL

    num_models <- length(models)

    for(i in 1:dim(scenarios)[1]){

      if(missing(arms)){
        arms <- c(2:scenarios[i,]$num_arms)
      }

      arms <- sort(arms)

      d_i <- as.numeric(scenarios[i, grepl("^d\\d", names(scenarios))])[1:scenarios[i,]$num_arms]
      theta_i <- as.numeric(scenarios[i, grepl("^theta\\d", names(scenarios))])[1:scenarios[i,]$num_arms]
      lambda_i <- as.numeric(scenarios[i, grepl("^lambda\\d", names(scenarios))])[1:(scenarios[i,]$num_arms+1)]

      time_dep_effect_i <- datasim_cont(n_arm = scenarios$n_arm[i],
                                        num_arms = scenarios$num_arms[i],
                                        d = d_i,
                                        period_blocks = scenarios$period_blocks[i],
                                        mu0 = scenarios$mu0[i],
                                        theta = theta_i,
                                        lambda = lambda_i,
                                        sigma = scenarios$sigma[i],
                                        trend = scenarios$trend[i],
                                        N_peak = scenarios$N_peak[i],
                                        n_wave = scenarios$n_wave[i],
                                        full = TRUE,
                                        check = FALSE)$time_dep_effect

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
                                                     n_wave = scenarios$n_wave[i],
                                                     full = FALSE,
                                                     check = FALSE),
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
                                 prec_theta = scenarios$prec_theta[i],
                                 prec_eta = scenarios$prec_eta[i],
                                 tau_a = scenarios$tau_a[i],
                                 tau_b = scenarios$tau_b[i],
                                 prec_a = scenarios$prec_a[i],
                                 prec_b = scenarios$prec_b[i],
                                 bucket_size = scenarios$bucket_size[i],
                                 smoothing_basis = scenarios$smoothing_basis[i],
                                 basis_dim = scenarios$basis_dim[i],
                                 gam_method = scenarios$gam_method[i],
                                 bs_degree = scenarios$bs_degree[i],
                                 poly_degree = scenarios$poly_degree[i])

                    }



      result_i <- cbind(scenarios[i,],
                        study_arm = rep(arms, each = num_models),
                        model = models,
                        reject_h0 = rowMeans(matrix(as.logical(unlist(unname(db[grep("reject_h0_", rownames(db)),]))), ncol = nsim), na.rm = TRUE), # get power/T1E
                        bias = rowMeans(matrix(as.double(unlist(unname(db[grep("treat_effect", rownames(db)),]))), ncol = nsim)-time_dep_effect_i[arms], na.rm = TRUE), # get bias (using time dependent treatment effect)
                        MSE = rowMeans((matrix(as.double(unlist(unname(db[grep("treat_effect", rownames(db)),]))), ncol = nsim)-time_dep_effect_i[arms])^2, na.rm = TRUE), # get MSE (using time dependent treatment effect)
                        failed = rowSums(is.na(matrix(as.logical(unlist(unname(db[grep("reject_h0_", rownames(db)),]))), ncol = nsim))),
                        nsim = nsim,
                        row.names = NULL)

      result <- rbind(result, result_i)

      print(paste0("Scenario ", i, "/", dim(scenarios)[1], " done. Time: ", Sys.time()))
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

      d_i <- as.numeric(scenarios[i, grepl("^d\\d", names(scenarios))])[1:scenarios[i,]$num_arms]
      OR_i <- as.numeric(scenarios[i, grepl("^OR\\d", names(scenarios))])[1:scenarios[i,]$num_arms]
      lambda_i <- as.numeric(scenarios[i, grepl("^lambda\\d", names(scenarios))])[1:(scenarios[i,]$num_arms+1)]

      time_dep_effect_i <- datasim_bin(n_arm = scenarios$n_arm[i],
                                       num_arms = scenarios$num_arms[i],
                                       d = d_i,
                                       period_blocks = scenarios$period_blocks[i],
                                       p0 = scenarios$p0[i],
                                       OR = OR_i,
                                       lambda = lambda_i,
                                       trend = scenarios$trend[i],
                                       N_peak = scenarios$N_peak[i],
                                       n_wave = scenarios$n_wave[i],
                                       full = TRUE,
                                       check = FALSE)$time_dep_effect

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
                                                    n_wave = scenarios$n_wave[i],
                                                    full = FALSE,
                                                    check = FALSE),
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
                                 prec_theta = scenarios$prec_theta[i],
                                 prec_eta = scenarios$prec_eta[i],
                                 tau_a = scenarios$tau_a[i],
                                 tau_b = scenarios$tau_b[i],
                                 prec_a = scenarios$prec_a[i],
                                 prec_b = scenarios$prec_b[i],
                                 bucket_size = scenarios$bucket_size[i],
                                 smoothing_basis = scenarios$smoothing_basis[i],
                                 basis_dim = scenarios$basis_dim[i],
                                 gam_method = scenarios$gam_method[i],
                                 bs_degree = scenarios$bs_degree[i],
                                 poly_degree = scenarios$poly_degree[i])
                    }



      result_i <- cbind(scenarios[i,],
                        study_arm = rep(arms, each = num_models),
                        model = models,
                        reject_h0 = rowMeans(matrix(as.logical(unlist(unname(db[grep("reject_h0_", rownames(db)),]))), ncol = nsim), na.rm = TRUE), # get power/T1E
                        bias = rowMeans(matrix(as.double(unlist(unname(db[grep("treat_effect", rownames(db)),]))), ncol = nsim)-log(time_dep_effect_i[arms]), na.rm = TRUE), # get bias (using time dependent treatment effect)
                        MSE = rowMeans((matrix(as.double(unlist(unname(db[grep("treat_effect", rownames(db)),]))), ncol = nsim)-log(time_dep_effect_i[arms]))^2, na.rm = TRUE), # get MSE (using time dependent treatment effect)
                        failed = rowSums(is.na(matrix(as.logical(unlist(unname(db[grep("reject_h0_", rownames(db)),]))), ncol = nsim))),
                        nsim = nsim,
                        row.names = NULL)

      result <- rbind(result, result_i)

      print(paste0("Scenario ", i, "/", dim(scenarios)[1], " done. Time: ", Sys.time()))

    }

    stopCluster(cl)
    gc()

  }
  return(result)
}



