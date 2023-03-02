#' Wrapper function performing simulation studies for a given set of scenarios (not parallelized)
#'
#' @description This function performs a simulation study for a given set of scenarios, analyzing simulated data using different models as indicated by the user. Performs inference for indicated experimental treatment arms. Simulates the probability to reject \eqn{H_0} based on a given number of replications.
#'
#' @param nsim Integer. Number of replications. Must be larger than 1.
#' @param scenarios Data frame containing all parameters for scenarios that should be simulated.
#' @param arms Integer vector with treatment arms to perform inference on. These arms are compared to the control group. Default - all arms except the first one.
#' @param models Character vector with models that should be used for the analysis. Default=c("fixmodel", "sepmodel", "poolmodel"). Available models for continuous endpoints are: 'fixmodel', 'fixmodel_cal', 'gam', 'MAPprior', 'mixmodel', 'mixmodel_cal', 'mixmodel_AR1', 'mixmodel_AR1_cal', 'piecewise', 'piecewise_cal', 'poolmodel', 'sepmodel', 'sepmodel_adj', 'splines', 'splines_cal', 'timemachine'. Available models for binary endpoints are: 'fixmodel', 'fixmodel_cal', 'MAPprior', 'poolmodel', 'sepmodel', 'sepmodel_adj', 'timemachine'.
#' @param endpoint Endpoint indicator. "cont" for continuous endpoints, "bin" for binary endpoints.
#' @param verbose Logical. Indicates whether to print a message (system time and number of finished scenarios) after simulating each scenario in order to track the progress of the simulations. Default=TRUE.
#'
#' @export
#'
#' @examples
#'
#' \donttest{
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
#' sim_results <- sim_study(nsim = 100, scenarios = sim_scenarios, arms = c(3, 4),
#' models = c("fixmodel", "sepmodel", "poolmodel"), endpoint = "cont")
#' }
#'
#'
#' @return Data frame with all considered scenarios and corresponding results - the probability to reject \eqn{H_0}.
#' @author Pavla Krotka


sim_study <- function(nsim, scenarios, arms, models = c("fixmodel", "sepmodel", "poolmodel"), endpoint, verbose = TRUE){

  if(!is.numeric(nsim) | length(nsim)!=1 | nsim<=1){
    stop("Number of replications (`nsim`) must be one number and must be larger than 1!")
  }

  if(!is.data.frame(scenarios)){
    stop("The considered scenarios (`scenarios`) must be a data frame!")
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



  # Add default parameters

  if(("period_blocks" %in% colnames(scenarios))==FALSE){
    scenarios$period_blocks <- 2
  }

  if(endpoint=="cont" & ("mu0" %in% colnames(scenarios))==FALSE){
    scenarios$mu0 <- 0
  }

  if(("alpha" %in% colnames(scenarios))==FALSE){
    scenarios$alpha <- 0.025
  }

  if((("fixmodel" %in% models) | ("fixmodel_cal" %in% models) | ("mixmodel" %in% models) | ("mixmodel_cal" %in% models) | ("mixmodel_AR1" %in% models) | ("mixmodel_AR1_cal" %in% models) | ("piecewise" %in% models) | ("piecewise_cal" %in% models) | ("splines" %in% models) | ("splines_cal" %in% models)) & ("ncc" %in% colnames(scenarios))==FALSE){
    scenarios$ncc <- TRUE
  }

  if((("fixmodel_cal" %in% models) | ("mixmodel_cal" %in% models) | ("mixmodel_AR1_cal" %in% models) | ("piecewise_cal" %in% models) | ("splines_cal" %in% models)) & ("unit_size" %in% colnames(scenarios))==FALSE){
    scenarios$unit_size <- 25
  }

  if((("mixmodel" %in% models) | ("mixmodel_cal" %in% models) | ("mixmodel_AR1" %in% models) | ("mixmodel_AR1_cal" %in% models) | ("gam" %in% models)) & ("ci" %in% colnames(scenarios))==FALSE){
    scenarios$ci <- FALSE
  }

  if(("gam" %in% models) & ("smoothing_basis" %in% colnames(scenarios))==FALSE){
    scenarios$smoothing_basis <- "tp"
  }

  if(("gam" %in% models) & ("basis_dim" %in% colnames(scenarios))==FALSE){
    scenarios$basis_dim <- -1
  }

  if(("gam" %in% models) & ("gam_method" %in% colnames(scenarios))==FALSE){
    scenarios$gam_method <- "GCV.Cp"
  }

  if((("splines" %in% models) | ("splines_cal" %in% models)) & ("bs_degree" %in% colnames(scenarios))==FALSE){
    scenarios$bs_degree <- 3
  }

  if((("piecewise" %in% models) | ("piecewise_cal" %in% models)) & ("poly_degree" %in% colnames(scenarios))==FALSE){
    scenarios$poly_degree <- 3
  }

  if(("MAPprior" %in% models) & ("opt" %in% colnames(scenarios))==FALSE){
    scenarios$opt <- 2
  }

  if(("MAPprior" %in% models) & ("prior_prec_tau" %in% colnames(scenarios))==FALSE){
    scenarios$prior_prec_tau <- 4
  }

  if(("MAPprior" %in% models) & ("prior_prec_eta" %in% colnames(scenarios))==FALSE){
    scenarios$prior_prec_eta <- 0.001
  }

  if(("MAPprior" %in% models) & ("n_samples" %in% colnames(scenarios))==FALSE){
    scenarios$n_samples <- 1000
  }

  if(("MAPprior" %in% models) & ("n_chains" %in% colnames(scenarios))==FALSE){
    scenarios$n_chains <- 4
  }

  if(("MAPprior" %in% models) & ("n_iter" %in% colnames(scenarios))==FALSE){
    scenarios$n_iter <- 4000
  }

  if(("MAPprior" %in% models) & ("n_adapt" %in% colnames(scenarios))==FALSE){
    scenarios$n_adapt <- 1000
  }

  if(("MAPprior" %in% models) & ("robustify" %in% colnames(scenarios))==FALSE){
    scenarios$robustify <- TRUE
  }

  if(("MAPprior" %in% models) & ("weight" %in% colnames(scenarios))==FALSE){
    scenarios$weight <- 0.1
  }

  if(("timemachine" %in% models) & ("prec_theta" %in% colnames(scenarios))==FALSE){
    scenarios$prec_theta <- 0.001
  }

  if(("timemachine" %in% models) & ("prec_eta" %in% colnames(scenarios))==FALSE){
    scenarios$prec_eta <- 0.001
  }

  if(("timemachine" %in% models) & ("tau_a" %in% colnames(scenarios))==FALSE){
    scenarios$tau_a <- 0.1
  }

  if(("timemachine" %in% models) & ("tau_b" %in% colnames(scenarios))==FALSE){
    scenarios$tau_b <- 0.01
  }

  if(("timemachine" %in% models) & ("bucket_size" %in% colnames(scenarios))==FALSE){
    scenarios$bucket_size <- 25
  }

  if(("timemachine" %in% models) & endpoint=="cont" & ("prec_a" %in% colnames(scenarios))==FALSE){
    scenarios$prec_a <- 0.001
  }

  if(("timemachine" %in% models) & endpoint=="cont" & ("prec_b" %in% colnames(scenarios))==FALSE){
    scenarios$prec_b <- 0.001
  }


  # Check scenarios data frame


  if(endpoint=="cont" & sum(c("num_arms", "n_arm", "d1", "theta1", "lambda1", "sigma", "trend") %in% colnames(scenarios)==FALSE)>0){
    stop("The `scenarios` data frame must include the parameters 'num_arms', 'n_arm', 'd', 'theta', 'lambda', 'sigma' and 'trend'!")
  }

  if(endpoint=="bin" & sum(c("num_arms", "n_arm", "d1", "p0", "OR1", "lambda1", "trend") %in% colnames(scenarios)==FALSE)>0){
    stop("The `scenarios` data frame must include the parameters 'num_arms', 'n_arm', 'd', 'p0', 'OR', 'lambda' and 'trend'!")
  }

  if(sum(scenarios$trend %in% c("linear", "linear_2", "stepwise", "stepwise_2", "inv_u", "seasonal")==FALSE)>0){
    stop("Values allowed for the time trend pattern (column 'trend') are: 'linear', 'linear_2', 'stepwise', 'stepwise_2', 'inv_u', 'seasonal'!")
  }

  if(("inv_u" %in% scenarios$trend) & ("N_peak" %in% colnames(scenarios))==FALSE){
    stop("If the time trend pattern is 'inv_u', the parameter 'N_peak' must be specified!")
  }

  if(("seasonal" %in% scenarios$trend) & ("n_wave" %in% colnames(scenarios))==FALSE){
    stop("If the time trend pattern is 'seasonal', the parameter 'n_wave' must be specified!")
  }

  if(sum(grepl("^d\\d", colnames(scenarios)))!=max(scenarios$num_arms)){
    stop("The number of columns specifying the parameter 'd' (columns must be named 'd1', 'd2', ect.) must correspond to the number of treatment arms ('n_arms')!")
  }

  if(endpoint=="cont" & sum(grepl("^theta\\d", colnames(scenarios)))!=max(scenarios$num_arms)){
    stop("The number of columns specifying the parameter 'theta' (columns must be named 'theta1', 'theta2', ect.) must correspond to the number of treatment arms ('n_arms')!")
  }

  if(endpoint=="bin" & sum(grepl("^OR\\d", colnames(scenarios)))!=max(scenarios$num_arms)){
    stop("The number of columns specifying the parameter 'OR' (columns must be named 'OR1', 'OR2', ect.) must correspond to the number of treatment arms ('n_arms')!")
  }

  if(sum(grepl("^lambda\\d", colnames(scenarios)))!=max(scenarios$num_arms)+1){
    stop("The number of columns specifying the parameter 'lambda' (columns must be named 'lambda0', 'lambda1', ect.) must correspond to the number of arms ('n_arms'+1)!")
  }





  if (verbose) {
    print(paste0("Starting the simulations. ", dim(scenarios)[1], " scenarios will be simulated. Starting time: ", Sys.time()))
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

      # d_i <- as.numeric(scenarios[i, grepl("^d\\d", names(scenarios))])

      d_i <- scenarios[i, grepl("^d\\d", names(scenarios))]

      if (sum(c(1:scenarios[i,]$num_arms) %in% as.numeric(gsub("d","", colnames(d_i))))!=scenarios[i,]$num_arms) {
        stop("Each experimental treatment arm needs to have an entry time specified in a separate column (columns must be named 'd1', 'd2', ect.)!")
      }

      d_i <- d_i[, order(as.numeric(gsub("d","", colnames(d_i))))] # Sort to correct order
      d_i <- as.numeric(d_i)[1:scenarios[i,]$num_arms]

      # theta_i <- as.numeric(scenarios[i, grepl("^theta\\d", names(scenarios))])

      theta_i <- scenarios[i, grepl("^theta\\d", names(scenarios))]

      if (sum(c(1:scenarios[i,]$num_arms) %in% as.numeric(gsub("theta","", colnames(theta_i))))!=scenarios[i,]$num_arms) {
        stop("Each experimental treatment arm needs to have a treatment effect specified in a separate column (columns must be named 'theta1', 'theta2', ect.)!")
      }

      theta_i <- theta_i[, order(as.numeric(gsub("theta","", colnames(theta_i))))] # Sort to correct order
      theta_i <- as.numeric(theta_i)[1:scenarios[i,]$num_arms]

      # lambda_i <- as.numeric(scenarios[i, grepl("^lambda\\d", names(scenarios))])

      lambda_i <- scenarios[i, grepl("^lambda\\d", names(scenarios))]

      if (sum(c(0:scenarios[i,]$num_arms) %in% as.numeric(gsub("lambda","", colnames(lambda_i))))!=(scenarios[i,]$num_arms+1)) {
        stop("Each arm needs to have a strength of time trend specified in a separate column (columns must be named 'lambda0', 'lambda1', ect.)!")
      }

      lambda_i <- lambda_i[, order(as.numeric(gsub("lambda","", colnames(lambda_i))))] # Sort to correct order
      lambda_i <- as.numeric(lambda_i)[1:(scenarios[i,]$num_arms+1)]

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
                                 prior_prec_eta = scenarios$prior_prec_eta[i],
                                 n_samples = scenarios$n_samples[i],
                                 n_chains = scenarios$n_chains[i],
                                 n_iter = scenarios$n_iter[i],
                                 n_adapt = scenarios$n_adapt[i],
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
                                 poly_degree = scenarios$poly_degree[i]))

      result_i <- cbind(scenarios[i,],
                        study_arm = rep(arms, each = num_models),
                        model = models,
                        reject_h0 = rowMeans(matrix(as.logical(unlist(unname(db[grep("reject_h0_", rownames(db)),]))), ncol = nsim), na.rm = TRUE), # get power/T1E
                        failed = rowSums(is.na(matrix(as.logical(unlist(unname(db[grep("reject_h0_", rownames(db)),]))), ncol = nsim))),
                        nsim = nsim,
                        row.names = NULL)

      result <- rbind(result, result_i)


      if (verbose) {
        print(paste0("Scenario ", i, "/", dim(scenarios)[1], " done. Time: ", Sys.time()))
      }

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

      # d_i <- as.numeric(scenarios[i, grepl("^d\\d", names(scenarios))])

      d_i <- scenarios[i, grepl("^d\\d", names(scenarios))]

      if (sum(c(1:scenarios[i,]$num_arms) %in% as.numeric(gsub("d","", colnames(d_i))))!=scenarios[i,]$num_arms) {
        stop("Each experimental treatment arm needs to have an entry time specified in a separate column (columns must be named 'd1', 'd2', ect.)!")
      }

      d_i <- d_i[, order(as.numeric(gsub("d","", colnames(d_i))))] # Sort to correct order
      d_i <- as.numeric(d_i)[1:scenarios[i,]$num_arms]

      # OR_i <- as.numeric(scenarios[i, grepl("^OR\\d", names(scenarios))])

      OR_i <- scenarios[i, grepl("^OR\\d", names(scenarios))]

      if (sum(c(1:scenarios[i,]$num_arms) %in% as.numeric(gsub("OR","", colnames(OR_i))))!=scenarios[i,]$num_arms) {
        stop("Each experimental treatment arm needs to have an odds ratio specified in a separate column (columns must be named 'OR1', 'OR2', ect.)!")
      }

      OR_i <- OR_i[, order(as.numeric(gsub("OR","", colnames(OR_i))))] # Sort to correct order
      OR_i <- as.numeric(OR_i)[1:scenarios[i,]$num_arms]

      # lambda_i <- as.numeric(scenarios[i, grepl("^lambda\\d", names(scenarios))])

      lambda_i <- scenarios[i, grepl("^lambda\\d", names(scenarios))]

      if (sum(c(0:scenarios[i,]$num_arms) %in% as.numeric(gsub("lambda","", colnames(lambda_i))))!=(scenarios[i,]$num_arms+1)) {
        stop("Each arm needs to have a strength of time trend specified in a separate column (columns must be named 'lambda0', 'lambda1', ect.)!")
      }

      lambda_i <- lambda_i[, order(as.numeric(gsub("lambda","", colnames(lambda_i))))] # Sort to correct order
      lambda_i <- as.numeric(lambda_i)[1:(scenarios[i,]$num_arms+1)]

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
                                 prior_prec_eta = scenarios$prior_prec_eta[i],
                                 n_samples = scenarios$n_samples[i],
                                 n_chains = scenarios$n_chains[i],
                                 n_iter = scenarios$n_iter[i],
                                 n_adapt = scenarios$n_adapt[i],
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
                                 poly_degree = scenarios$poly_degree[i]))

      result_i <- cbind(scenarios[i,],
                        study_arm = rep(arms, each = num_models),
                        model = models,
                        reject_h0 = rowMeans(matrix(as.logical(unlist(unname(db[grep("reject_h0_", rownames(db)),]))), ncol = nsim), na.rm = TRUE), # get power/T1E
                        failed = rowSums(is.na(matrix(as.logical(unlist(unname(db[grep("reject_h0_", rownames(db)),]))), ncol = nsim))),
                        nsim = nsim,
                        row.names = NULL)

      result <- rbind(result, result_i)


      if (verbose) {
        print(paste0("Scenario ", i, "/", dim(scenarios)[1], " done. Time: ", Sys.time()))
      }

    }
  }
  return(result)
}



